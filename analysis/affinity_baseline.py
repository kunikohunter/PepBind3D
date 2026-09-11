"""
Reviewer 1, major point 3 (RESPONSE_TO_REVIEWERS.md): "develop a simple
ML-based affinity predictor based on the near-native structure ensembles
and/or Rosetta scores to validate the value of PepBind3D". Implements the
pre-registered specification at
/data/p_csb_meiler/huntek1/revision/BASELINE_MODEL_SPEC.md (written before
any number was computed; section 8 fixes the reporting language in advance).

This is a Technical Validation check, not a modelling contribution: a
deliberately minimal, fixed-hyperparameter gradient-boosted-trees baseline on
the tabular Rosetta score columns already shipped in metadata.csv, evaluated
under leakage-controlled, grouped cross-validation, against a sequence-only
reference trained on identical folds. Per spec section 2.4, this baseline is
architecturally restricted to tabular scores: no graph/equivariant/learned
structural representation is used anywhere in this file, and none should be
added -- a separate, unpublished E(n)-equivariant structure-learning study is
in preparation by the same author, and this descriptor must not preempt it.

Censoring constants and the censoring rule are imported (not redefined) from
censored_vs_quantitative_auroc.py so "censored" means one thing across the
paper (spec section 1.1, pre-flight checklist item 1). IC50 and Kd are always
modelled as two separate models (spec 1.1); Kd reporting is additionally
stratified by the recovered raw IEDB sub-label per spec 1.4, using the join
already computed by kd_label_pooling.py (kd_label_recovery.csv) -- no
retraining.

Usage:
    $ python3 affinity_baseline.py --out-dir <output dir>

Self-test (synthetic-data recovery + the CV-1 vs naive-random-split leakage
test that proves the grouped-split machinery actually removes leakage):
    $ python3 affinity_baseline.py --self-test
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score, mean_squared_error, roc_auc_score
from sklearn.model_selection import GroupKFold, KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from censored_vs_quantitative_auroc import (  # noqa: E402
    HF_DIR,
    IC50_CEILINGS,
    IC50_FLOOR,
    KD_CEILINGS,
    KD_FLOOR,
    STRONG_BINDER_THRESHOLD_NM,
    is_censored,
)

METADATA_FN = HF_DIR / "metadata.csv"
KD_LABEL_RECOVERY_FN = Path(__file__).resolve().parent / "results" / "kd_label_recovery.csv"

COL_ALLELE = "allele"
COL_ALLELE_IEDB = "allele_iedb"
COL_PEPTIDE = "peptide"
COL_PEPLEN = "peptide_length"
COL_MEAS_TYPE = "measurement_type"
COL_MEAS_VALUE = "measurement_value"
COL_PUBMED = "pubmed_id"
COL_FLAGGED = "flagged"

# --- spec section 2.1: the 9 structural (Rosetta) features -----------------
STRUCT_SHIPPED = [
    "I_sc_best", "I_sc_mean",
    "reweighted_sc_best", "reweighted_sc_mean",
    "total_score_best", "total_score_mean",
]
STRUCT_GAP = ["I_sc_gap", "reweighted_sc_gap", "total_score_gap"]
STRUCT_FEATURES = STRUCT_SHIPPED + STRUCT_GAP  # 9 columns, exactly per spec

# rosetta_best_score/rosetta_mean_score are exact duplicates of
# total_score_best/total_score_mean in the released file (spec 2.1); dropped
# after asserting the duplication, never used as features.
DUPLICATE_COLS = {"rosetta_best_score": "total_score_best", "rosetta_mean_score": "total_score_mean"}

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
N_SEQ_SLOTS = 9  # spec 5: first-4 + last-5 anchor-preserving alignment

ASSAYS = [
    ("IC50", "IC50", IC50_CEILINGS, IC50_FLOOR),
    ("Kd", "Kd", KD_CEILINGS, KD_FLOOR),
]

# Fixed hyperparameters, spec section 3. Not tuned. No search, no ensembling,
# no seed averaging, no calibration, no feature selection.
HGB_KWARGS = dict(
    learning_rate=0.05,
    max_iter=300,
    max_leaf_nodes=31,
    min_samples_leaf=50,
    l2_regularization=1.0,
    early_stopping=False,
    random_state=0,
)

MIN_PER_ALLELE_N = 10        # spec 6: per-allele reporting threshold
QUOTED_PER_ALLELE_N = 100    # spec 6: alleles quoted in prose
MIN_ALLELE_FOLD_N = 100      # spec 4.3: min pairs of the relevant assay for a CV-2 fold
N_BOOT_POOLED = 1000
N_BOOT_PER_ALLELE = 300


# =============================================================================
# Union-find (peptide clustering, spec 4.2)
# =============================================================================
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _hamming_link_same_length(peptides_by_len, uf, index_of, chunk=1500):
    """Rule (a): same length, Hamming distance <= 2. Vectorized, chunked over
    rows so memory stays bounded even for the ~10k-peptide 9-mer block."""
    for length, peps in peptides_by_len.items():
        n = len(peps)
        if n < 2:
            continue
        arr = np.array([[ord(c) for c in p] for p in peps], dtype=np.uint8)  # (n, length)
        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            block = arr[start:end]  # (b, length)
            # mismatches[i, j] = number of differing positions between
            # block[i] and arr[j], for j >= start (upper triangle only,
            # avoids double work and self-comparison)
            mism = np.zeros((end - start, n - start), dtype=np.int16)
            for pos in range(length):
                mism += (block[:, pos][:, None] != arr[start:, pos][None, :])
            close = np.argwhere((mism <= 2))
            for i, j in close:
                gi, gj = start + i, start + j
                if gj <= gi:
                    continue
                uf.union(index_of[peps[gi]], index_of[peps[gj]])


def _substring_link_cross_length(peptides_by_len, uf, index_of):
    """Rule (b): different length, shorter is a contiguous substring of the
    longer (N/C-terminal extension variants). Implemented as substring
    generation + hash-set membership, not all-pairs comparison, so it stays
    cheap regardless of the length gap."""
    lengths = sorted(peptides_by_len)
    length_sets = {length: set(peps) for length, peps in peptides_by_len.items()}
    for long_len in lengths:
        for long_pep in peptides_by_len[long_len]:
            for short_len in lengths:
                if short_len >= long_len:
                    continue
                short_set = length_sets[short_len]
                for i in range(long_len - short_len + 1):
                    sub = long_pep[i:i + short_len]
                    if sub in short_set:
                        uf.union(index_of[long_pep], index_of[sub])


def cluster_peptides(peptides):
    """Single-linkage clustering of unique peptide strings per spec 4.2.
    Returns (peptide -> cluster_id dict, diagnostics dict)."""
    uniq = sorted(set(peptides))
    index_of = {p: i for i, p in enumerate(uniq)}
    uf = UnionFind(len(uniq))

    by_len = {}
    for p in uniq:
        by_len.setdefault(len(p), []).append(p)

    _hamming_link_same_length(by_len, uf, index_of)
    _substring_link_cross_length(by_len, uf, index_of)

    roots = [uf.find(i) for i in range(len(uniq))]
    root_to_cluster = {}
    cluster_id = np.empty(len(uniq), dtype=np.int64)
    for i, r in enumerate(roots):
        if r not in root_to_cluster:
            root_to_cluster[r] = len(root_to_cluster)
        cluster_id[i] = root_to_cluster[r]

    pep_to_cluster = dict(zip(uniq, cluster_id))
    sizes = pd.Series(cluster_id).value_counts()
    diagnostics = {
        "n_unique_peptides": len(uniq),
        "n_clusters": int(len(sizes)),
        "largest_cluster_size": int(sizes.max()),
        "frac_singleton_clusters": float((sizes == 1).sum() / len(sizes)),
    }
    return pep_to_cluster, diagnostics


# =============================================================================
# Sequence one-hot features (spec section 5, B2/B3)
# =============================================================================
def _align9(peptide):
    """First-4 + last-5 anchor-preserving alignment into 9 slots (spec 5).
    Not a learned embedding; a fixed, deterministic re-indexing."""
    if len(peptide) <= N_SEQ_SLOTS:
        return (peptide[:4] + peptide[-5:])[:N_SEQ_SLOTS] if len(peptide) >= 5 else peptide.ljust(N_SEQ_SLOTS, "X")
    return peptide[:4] + peptide[-5:]


def seq_onehot(peptides):
    aa_index = {aa: i for i, aa in enumerate(AA20)}
    n = len(peptides)
    X = np.zeros((n, N_SEQ_SLOTS * len(AA20)), dtype=np.float32)
    for row, pep in enumerate(peptides):
        aligned = _align9(pep)
        for slot, ch in enumerate(aligned):
            j = aa_index.get(ch)
            if j is not None:
                X[row, slot * len(AA20) + j] = 1.0
    return X


# =============================================================================
# Data loading / filtering (spec 1.1, 2.1)
# =============================================================================
def load_and_prepare(metadata_fn=METADATA_FN):
    df = pd.read_csv(metadata_fn, low_memory=False)

    # pre-flight checklist item: assert (don't assume) the duplicate columns,
    # then drop them.
    for dup_col, canon_col in DUPLICATE_COLS.items():
        mismatch = (df[dup_col] - df[canon_col]).abs() > 1e-6
        mismatch &= ~(df[dup_col].isna() & df[canon_col].isna())
        assert mismatch.sum() == 0, (
            f"{dup_col} is not an exact duplicate of {canon_col} in {mismatch.sum()} rows; "
            "spec 2.1's drop-as-duplicate decision no longer holds, stop and re-check."
        )
    df = df.drop(columns=list(DUPLICATE_COLS.keys()))

    df = df[~df[COL_FLAGGED].astype(bool)].copy()

    df["I_sc_gap"] = df["I_sc_mean"] - df["I_sc_best"]
    df["reweighted_sc_gap"] = df["reweighted_sc_mean"] - df["reweighted_sc_best"]
    df["total_score_gap"] = df["total_score_mean"] - df["total_score_best"]

    return df


def filter_assay(df, meas_type, ceilings, floor):
    """Rows for one assay type: flagged already dropped upstream; non-null
    measurement_value and I_sc_best; not censored (spec 1.1)."""
    sub = df[df[COL_MEAS_TYPE] == meas_type].copy()
    sub = sub.dropna(subset=[COL_MEAS_VALUE, "I_sc_best"])
    censored = is_censored(sub[COL_MEAS_VALUE], ceilings, floor)
    sub = sub[~censored].copy()

    # Spec 4.4 asserts curation already deduplicates per (allele, peptide) and
    # says the script should assert this rather than assume it. On the
    # released file this holds for Kd but NOT for IC50: one (allele, peptide)
    # pair (A*02:01 / MLYQLLEAV, PMID 19734234) has two IC50 rows with the
    # same Rosetta scores (same threaded structure) but different reported
    # measurement_value (14949 nM vs 8.2 nM) -- an un-deduplicated curation
    # artefact, not a script bug. Rather than crash, keep the first
    # occurrence (deterministic) and record the finding instead of silently
    # assuming the spec's premise. This affects 1 of ~10,900 IC50 rows and is
    # reported in affinity_baseline_metrics.json.
    dup_mask = sub.duplicated(subset=[COL_ALLELE, COL_PEPTIDE], keep=False)
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows:
        print(f"  WARNING: {n_dup_rows} rows involved in (allele, peptide) duplicates for {meas_type} "
              "(spec 4.4 assumed curation deduplicates; keeping first occurrence per pair).")
    sub = sub.drop_duplicates(subset=[COL_ALLELE, COL_PEPTIDE], keep="first")
    sub.attrs["n_duplicate_rows_dropped"] = n_dup_rows

    sub["log10_affinity"] = np.log10(sub[COL_MEAS_VALUE])
    sub["strong_binder"] = (sub[COL_MEAS_VALUE] < STRONG_BINDER_THRESHOLD_NM).astype(int)
    sub[COL_ALLELE] = sub[COL_ALLELE].astype("category")
    sub["allele_code"] = sub[COL_ALLELE].cat.codes.astype(np.float64)
    return sub.reset_index(drop=True)


# =============================================================================
# Feature assembly
# =============================================================================
def build_features(sub, include_allele):
    """Returns (X_struct_only, X_seq_only, X_full, allele_col_index_in_full)
    as numpy arrays; allele is appended last when included."""
    struct = sub[STRUCT_FEATURES].to_numpy(dtype=np.float64)
    seq = seq_onehot(sub[COL_PEPTIDE].tolist())
    peplen = sub[[COL_PEPLEN]].to_numpy(dtype=np.float64)

    seq_block = np.hstack([seq, peplen])
    if include_allele:
        allele_col = sub[["allele_code"]].to_numpy(dtype=np.float64)
        # B1: unfitted single score. I_sc_best is used directly (not negated) --
        # verified empirically (see report) that raw I_sc_best is POSITIVELY
        # correlated with log10(affinity) in this dataset (both lower REU and
        # lower affinity/nM track tighter binding), matching the manuscript's
        # Figure 3 convention. An earlier draft negated I_sc_best by analogy with
        # censored_vs_quantitative_auroc.py's AUROC scoring convention, which does
        # not apply here and silently flipped the sign of every B1 metric and every
        # B3-B1 delta -- caught by cross-checking against a raw correlation before
        # accepting the first full run's numbers.
        X_b1 = sub[["I_sc_best"]].to_numpy(dtype=np.float64)
        X_b2 = np.hstack([seq_block, allele_col])
        X_b3 = np.hstack([seq_block, struct, allele_col])
        allele_idx_b2 = X_b2.shape[1] - 1
        allele_idx_b3 = X_b3.shape[1] - 1
    else:
        X_b1 = sub[["I_sc_best"]].to_numpy(dtype=np.float64)
        X_b2 = seq_block
        X_b3 = np.hstack([seq_block, struct])
        allele_idx_b2 = None
        allele_idx_b3 = None
    return {"B0": X_b1, "B1": X_b1, "B2": X_b2, "B3": X_b3,
            "allele_idx": {"B2": allele_idx_b2, "B3": allele_idx_b3}}


# =============================================================================
# Split schemes (spec section 4)
# =============================================================================
def cv1_splits(sub, n_splits=5):
    """GroupKFold, groups = peptide cluster. allele stays in the feature
    set (shared between train/test)."""
    groups = sub["cluster_id"].to_numpy()
    n_groups = len(np.unique(groups))
    gkf = GroupKFold(n_splits=min(n_splits, n_groups))
    for train_idx, test_idx in gkf.split(sub, groups=groups):
        yield train_idx, test_idx, None


def cv2_splits(sub, min_n=MIN_ALLELE_FOLD_N):
    """Leave-one-allele-out. Only alleles with >= min_n rows become a test
    fold; smaller alleles stay in every training fold and are never tested
    (spec 4.3). allele is NOT a usable feature here (dropped by the caller)."""
    counts = sub[COL_ALLELE].value_counts()
    qualifying = counts[counts >= min_n].index.tolist()
    idx_by_allele = {a: sub.index[sub[COL_ALLELE] == a].to_numpy() for a in sub[COL_ALLELE].unique()}
    all_idx = sub.index.to_numpy()
    for allele in qualifying:
        test_idx = idx_by_allele[allele]
        train_idx = np.setdiff1d(all_idx, test_idx, assume_unique=False)
        yield train_idx, test_idx, allele


def cv3_splits(sub, n_blocks=5, seed=0):
    """Optional strictest scheme: allele blocks x peptide-cluster blocks
    (spec 4.3). Rows in exactly one held-out block are discarded for that
    fold. allele is not usable (held-out allele-block rows never enter train
    for that fold -- see affinity_baseline docstring reasoning)."""
    rng = np.random.default_rng(seed)
    alleles = sub[COL_ALLELE].astype(str).unique()
    clusters = sub["cluster_id"].unique()
    rng.shuffle(alleles)
    rng.shuffle(clusters)
    allele_block = {a: i % n_blocks for i, a in enumerate(alleles)}
    cluster_block = {c: i % n_blocks for i, c in enumerate(clusters)}

    allele_block_arr = sub[COL_ALLELE].astype(str).map(allele_block).to_numpy()
    cluster_block_arr = sub["cluster_id"].map(cluster_block).to_numpy()

    for k in range(n_blocks):
        in_allele_k = allele_block_arr == k
        in_cluster_k = cluster_block_arr == k
        test_mask = in_allele_k & in_cluster_k
        train_mask = (~in_allele_k) & (~in_cluster_k)
        test_idx = sub.index[test_mask].to_numpy()
        train_idx = sub.index[train_mask].to_numpy()
        if len(test_idx) == 0 or len(train_idx) == 0:
            continue
        yield train_idx, test_idx, None


# =============================================================================
# Models B0-B3
# =============================================================================
def fit_predict_regression(model_name, X_train, y_train, X_test, allele_idx):
    if model_name == "B0":
        pred = np.full(len(X_test), float(np.mean(y_train)))
        return pred
    if model_name == "B1":
        # X_test here is already -I_sc_best; unfitted, ranking-only feature.
        return X_test[:, 0]
    categorical_features = [allele_idx] if allele_idx is not None else None
    reg = HistGradientBoostingRegressor(categorical_features=categorical_features, **HGB_KWARGS)
    reg.fit(X_train, y_train)
    return reg.predict(X_test)


def fit_predict_classification(model_name, X_train, y_train, X_test, allele_idx):
    if model_name == "B0":
        prevalence = float(np.mean(y_train))
        return np.full(len(X_test), prevalence)
    if model_name == "B1":
        return X_test[:, 0]
    categorical_features = [allele_idx] if allele_idx is not None else None
    if len(np.unique(y_train)) < 2:
        return np.full(len(X_test), float(np.mean(y_train)))
    clf = HistGradientBoostingClassifier(categorical_features=categorical_features, **HGB_KWARGS)
    clf.fit(X_train, y_train)
    return clf.predict_proba(X_test)[:, 1]


# =============================================================================
# Metrics + cluster-level bootstrap (spec section 6)
# =============================================================================
def _safe_spearman(y_true, y_pred):
    if len(y_true) < 3 or np.std(y_pred) == 0:
        return np.nan
    return spearmanr(y_true, y_pred).correlation


def _safe_pearson(y_true, y_pred):
    if len(y_true) < 3 or np.std(y_pred) == 0:
        return np.nan
    return pearsonr(y_true, y_pred)[0]


def _rmse(y_true, y_pred):
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


def bootstrap_ci_clustered(cluster_ids, metric_fn, n_boot=N_BOOT_POOLED, seed=0):
    """Percentile bootstrap CI, resampling at the peptide-cluster level (spec
    6): a row bootstrap on grouped data understates the true CI. Mirrors the
    skip-invalid-resample / require-half-to-survive idiom of
    censored_vs_quantitative_auroc.bootstrap_auroc_ci."""
    rng = np.random.default_rng(seed)
    uniq_clusters = np.unique(cluster_ids)
    rows_by_cluster = {c: np.where(cluster_ids == c)[0] for c in uniq_clusters}
    n_clusters = len(uniq_clusters)
    boots = []
    for _ in range(n_boot):
        sampled = rng.choice(uniq_clusters, n_clusters, replace=True)
        idx = np.concatenate([rows_by_cluster[c] for c in sampled])
        val = metric_fn(idx)
        if val is not None and np.isfinite(val):
            boots.append(val)
    if len(boots) < n_boot // 2:
        return None, None
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def regression_metrics_with_ci(y_true, y_pred, cluster_ids, n_boot=N_BOOT_POOLED, seed=0):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    rho = _safe_spearman(y_true, y_pred)
    r = _safe_pearson(y_true, y_pred)
    rmse = _rmse(y_true, y_pred)

    def m_rho(idx):
        return _safe_spearman(y_true[idx], y_pred[idx])

    def m_r(idx):
        return _safe_pearson(y_true[idx], y_pred[idx])

    def m_rmse(idx):
        return _rmse(y_true[idx], y_pred[idx])

    rho_lo, rho_hi = bootstrap_ci_clustered(cluster_ids, m_rho, n_boot, seed)
    r_lo, r_hi = bootstrap_ci_clustered(cluster_ids, m_r, n_boot, seed + 1)
    rmse_lo, rmse_hi = bootstrap_ci_clustered(cluster_ids, m_rmse, n_boot, seed + 2)
    return {
        "n": int(len(y_true)),
        "spearman": None if pd.isna(rho) else float(rho), "spearman_ci_lo": rho_lo, "spearman_ci_hi": rho_hi,
        "pearson": None if pd.isna(r) else float(r), "pearson_ci_lo": r_lo, "pearson_ci_hi": r_hi,
        "rmse": rmse, "rmse_ci_lo": rmse_lo, "rmse_ci_hi": rmse_hi,
    }


def classification_metrics_with_ci(y_true, y_score, cluster_ids, n_boot=N_BOOT_POOLED, seed=0):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    prevalence = float(np.mean(y_true))
    if len(np.unique(y_true)) < 2:
        return {"n": int(len(y_true)), "prevalence": prevalence,
                "auroc": None, "auroc_ci_lo": None, "auroc_ci_hi": None,
                "auprc": None, "auprc_ci_lo": None, "auprc_ci_hi": None}
    auroc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))

    def m_auroc(idx):
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            return None
        return roc_auc_score(yt, y_score[idx])

    def m_auprc(idx):
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            return None
        return average_precision_score(yt, y_score[idx])

    auroc_lo, auroc_hi = bootstrap_ci_clustered(cluster_ids, m_auroc, n_boot, seed)
    auprc_lo, auprc_hi = bootstrap_ci_clustered(cluster_ids, m_auprc, n_boot, seed + 1)
    return {
        "n": int(len(y_true)), "prevalence": prevalence,
        "auroc": auroc, "auroc_ci_lo": auroc_lo, "auroc_ci_hi": auroc_hi,
        "auprc": auprc, "auprc_ci_lo": auprc_lo, "auprc_ci_hi": auprc_hi,
    }


def paired_delta_ci(y_true, pred_a, pred_b, cluster_ids, metric, n_boot=N_BOOT_POOLED, seed=0):
    """Bootstrap CI of the difference metric(pred_a) - metric(pred_b) on the
    SAME resampled clusters each iteration (spec 6, 'paired deltas')."""
    y_true = np.asarray(y_true)
    pred_a = np.asarray(pred_a)
    pred_b = np.asarray(pred_b)

    if metric == "spearman":
        point = _safe_spearman(y_true, pred_a) - _safe_spearman(y_true, pred_b)
    elif metric == "auroc":
        point = roc_auc_score(y_true, pred_a) - roc_auc_score(y_true, pred_b) if len(np.unique(y_true)) > 1 else np.nan
    else:
        raise ValueError(metric)

    def m(idx):
        yt = y_true[idx]
        if metric == "spearman":
            a = _safe_spearman(yt, pred_a[idx])
            b = _safe_spearman(yt, pred_b[idx])
            if a is None or b is None or np.isnan(a) or np.isnan(b):
                return None
            return a - b
        if len(np.unique(yt)) < 2:
            return None
        return roc_auc_score(yt, pred_a[idx]) - roc_auc_score(yt, pred_b[idx])

    lo, hi = bootstrap_ci_clustered(cluster_ids, m, n_boot, seed)
    return {"delta": None if pd.isna(point) else float(point), "ci_lo": lo, "ci_hi": hi}


def summarize_per_allele(rows, y_col, pred_col, allele_col, cluster_col, task):
    """Per-allele metrics with n>=10 threshold and median/IQR summary,
    matching summarize_per_allele() convention in
    02_score_affinity_validation.ipynb (spec 6)."""
    out = []
    for allele, g in rows.groupby(allele_col):
        if len(g) < MIN_PER_ALLELE_N:
            continue
        if task == "regression":
            m = regression_metrics_with_ci(g[y_col], g[pred_col], g[cluster_col].to_numpy(),
                                            n_boot=N_BOOT_PER_ALLELE, seed=0)
            out.append({"allele": allele, **m})
        else:
            m = classification_metrics_with_ci(g[y_col], g[pred_col], g[cluster_col].to_numpy(),
                                                n_boot=N_BOOT_PER_ALLELE, seed=0)
            out.append({"allele": allele, **m})
    per_allele_df = pd.DataFrame(out)
    summary = {}
    if task == "regression" and len(per_allele_df):
        rhos = per_allele_df["spearman"].dropna()
        summary = {
            "median_spearman": float(rhos.median()) if len(rhos) else None,
            "iqr_lo": float(rhos.quantile(0.25)) if len(rhos) else None,
            "iqr_hi": float(rhos.quantile(0.75)) if len(rhos) else None,
            "n_alleles": int(len(rhos)),
            "n_alleles_rho_gt_0.3": int((rhos > 0.3).sum()),
            "n_alleles_rho_lt_0": int((rhos < 0).sum()),
        }
    return per_allele_df, summary


# =============================================================================
# Orchestration
# =============================================================================
def run_scheme(sub, scheme_name, include_allele, task):
    """Fits B0-B3 across all folds of one scheme, returns concatenated
    out-of-fold predictions as a DataFrame (one row per held-out pair)."""
    feats = build_features(sub, include_allele=include_allele)
    y_reg = sub["log10_affinity"].to_numpy()
    y_clf = sub["strong_binder"].to_numpy()

    if scheme_name == "cv1":
        splitter = cv1_splits(sub)
    elif scheme_name == "cv2":
        splitter = cv2_splits(sub)
    elif scheme_name == "cv3":
        splitter = cv3_splits(sub)
    else:
        raise ValueError(scheme_name)

    records = []
    n_folds = 0
    for train_idx, test_idx, held_out_allele in splitter:
        n_folds += 1
        fold_pred = {}
        for model_name in ("B0", "B1", "B2", "B3"):
            X = feats[model_name]
            allele_idx = feats["allele_idx"].get(model_name)
            if task == "regression":
                pred = fit_predict_regression(model_name, X[train_idx], y_reg[train_idx], X[test_idx], allele_idx)
            else:
                pred = fit_predict_classification(model_name, X[train_idx], y_clf[train_idx], X[test_idx], allele_idx)
            fold_pred[model_name] = pred
        fold_df = pd.DataFrame({
            "fold": n_folds, "held_out_allele": held_out_allele,
            "allele": sub[COL_ALLELE].to_numpy()[test_idx],
            "cluster_id": sub["cluster_id"].to_numpy()[test_idx],
            "peptide": sub[COL_PEPTIDE].to_numpy()[test_idx],
            "allele_iedb": sub[COL_ALLELE_IEDB].to_numpy()[test_idx],
            "measurement_value": sub[COL_MEAS_VALUE].to_numpy()[test_idx],
            "pubmed_id": sub[COL_PUBMED].to_numpy()[test_idx],
            "y_reg": y_reg[test_idx], "y_clf": y_clf[test_idx],
        })
        for model_name in ("B0", "B1", "B2", "B3"):
            fold_df[f"pred_{model_name}"] = fold_pred[model_name]
        records.append(fold_df)
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame(), n_folds


def summarize_scheme(oof, task, assay_label, scheme_name):
    y_col = "y_reg" if task == "regression" else "y_clf"
    out = {"assay": assay_label, "scheme": scheme_name, "task": task, "n_test_rows": int(len(oof))}
    pooled = {}
    per_allele_all = {}
    for model_name in ("B0", "B1", "B2", "B3"):
        pred_col = f"pred_{model_name}"
        if task == "regression":
            pooled[model_name] = regression_metrics_with_ci(
                oof[y_col], oof[pred_col], oof["cluster_id"].to_numpy(), seed=hash(model_name) % 1000)
        else:
            pooled[model_name] = classification_metrics_with_ci(
                oof[y_col], oof[pred_col], oof["cluster_id"].to_numpy(), seed=hash(model_name) % 1000)
        per_allele_df, summ = summarize_per_allele(oof, y_col, pred_col, "allele", "cluster_id", task)
        per_allele_all[model_name] = (per_allele_df, summ)
    out["pooled"] = pooled
    out["per_allele_summary"] = {m: s for m, (_, s) in per_allele_all.items()}

    metric = "spearman" if task == "regression" else "auroc"
    deltas = {}
    if task == "regression":
        deltas["B3_minus_B2"] = paired_delta_ci(oof[y_col], oof["pred_B3"], oof["pred_B2"],
                                                 oof["cluster_id"].to_numpy(), metric, seed=11)
        deltas["B3_minus_B1"] = paired_delta_ci(oof[y_col], oof["pred_B3"], oof["pred_B1"],
                                                 oof["cluster_id"].to_numpy(), metric, seed=12)
    else:
        deltas["B3_minus_B2"] = paired_delta_ci(oof[y_col], oof["pred_B3"], oof["pred_B2"],
                                                 oof["cluster_id"].to_numpy(), metric, seed=13)
        deltas["B3_minus_B1"] = paired_delta_ci(oof[y_col], oof["pred_B3"], oof["pred_B1"],
                                                 oof["cluster_id"].to_numpy(), metric, seed=14)
    out["paired_deltas"] = deltas
    return out, {m: per_allele_all[m][0] for m in per_allele_all}


def kd_sublabel_sensitivity(oof_cv1_kd_b3):
    """Spec 1.4: pooling was tested and found not-justified by
    kd_label_pooling.py (the two present raw sub-labels differ, Kruskal-
    Wallis / KS p << 0.05; the third labeled variant never matches the
    curation pipeline's normalization regex and has n=0). No retraining --
    just stratify the already-computed CV-1 B3 Kd predictions by recovered
    raw label, and flag if one sub-label's metric differs materially."""
    if not KD_LABEL_RECOVERY_FN.exists():
        return {"note": f"{KD_LABEL_RECOVERY_FN} not found; skipping sub-label sensitivity."}
    recovery = pd.read_csv(KD_LABEL_RECOVERY_FN)
    for d in (oof_cv1_kd_b3, recovery):
        d["pubmed_id"] = d["pubmed_id"].astype(str).str.strip()
        d["peptide"] = d["peptide"].astype(str).str.strip()
        d["allele_iedb"] = d["allele_iedb"].astype(str).str.strip()
    recovery = recovery.drop_duplicates(subset=["allele_iedb", "peptide", "measurement_value", "pubmed_id"])
    merged = oof_cv1_kd_b3.merge(
        recovery[["allele_iedb", "peptide", "measurement_value", "pubmed_id", "raw_label"]],
        on=["allele_iedb", "peptide", "measurement_value", "pubmed_id"], how="left",
    )
    out = {}
    for label, g in merged.groupby("raw_label"):
        if len(g) < MIN_PER_ALLELE_N:
            continue
        m = regression_metrics_with_ci(g["y_reg"], g["pred_B3"], g["cluster_id"].to_numpy(),
                                        n_boot=N_BOOT_PER_ALLELE, seed=0)
        out[label] = {**m, "n": int(len(g))}
    out["n_matched"] = int(merged["raw_label"].notna().sum())
    out["n_total"] = int(len(merged))
    return out


def cluster_diagnostics_to_dict(diag, df_all):
    peptide_alleles = df_all.groupby(COL_PEPTIDE)[COL_ALLELE].nunique()
    multi_allele_peptides = (peptide_alleles > 1).sum()
    rows_with_multi_allele_peptide = df_all[COL_PEPTIDE].map(peptide_alleles > 1).sum()
    diag = dict(diag)
    diag["frac_peptides_multi_allele"] = float(multi_allele_peptides / len(peptide_alleles))
    diag["frac_rows_peptide_seen_with_other_allele"] = float(rows_with_multi_allele_peptide / len(df_all))
    return diag


# =============================================================================
# Self-test
# =============================================================================
def _random_peptides(rng, n, length=9):
    return ["".join(rng.choice(AA20, length)) for _ in range(n)]


def self_test():
    rng = np.random.default_rng(0)

    # (i) known analytic recovery: y = -0.9*x + noise, HGB regressor on x
    # alone should recover rho close to the analytic rank correlation of the
    # generating relationship.
    n = 2000
    x = rng.normal(size=n)
    noise = rng.normal(scale=0.15, size=n)
    y = -0.9 * x + noise
    expected_rho = spearmanr(x, y).correlation  # analytic target for this exact draw
    X = x.reshape(-1, 1)
    reg = HistGradientBoostingRegressor(**HGB_KWARGS)
    reg.fit(X[:1500], y[:1500])
    pred = reg.predict(X[1500:])
    rho = spearmanr(y[1500:], pred).correlation
    assert abs(rho - abs(expected_rho)) < 0.15 or rho > 0.8, (
        f"self-test FAILED: recovered rho {rho:.3f} not close to expected {expected_rho:.3f}"
    )

    # (ii) leakage test: label is a PURE function of peptide identity (no
    # signal anywhere else). A model given peptide identity as a feature must
    # score ~0 under CV-1 (peptide-cluster-grouped) and ~1 under a naive
    # random row split. This is the test that proves the split code works.
    n_peptides = 200
    replicates = 20
    peptides = _random_peptides(rng, n_peptides, length=9)
    pep_label = {p: int(rng.integers(0, 2)) for p in peptides}
    rows_pep = np.repeat(peptides, replicates)
    rng.shuffle(rows_pep)
    y_leak = np.array([pep_label[p] for p in rows_pep])
    pep_to_code = {p: i for i, p in enumerate(peptides)}
    peptide_code = np.array([pep_to_code[p] for p in rows_pep], dtype=np.float64).reshape(-1, 1)

    pep_to_cluster, _ = cluster_peptides(list(peptides))
    cluster_ids = np.array([pep_to_cluster[p] for p in rows_pep])

    # naive random split (row-level, ignores peptide identity)
    kf = KFold(n_splits=5, shuffle=True, random_state=0)
    naive_aucs = []
    for tr, te in kf.split(peptide_code):
        clf = HistGradientBoostingClassifier(categorical_features=[0], **HGB_KWARGS)
        clf.fit(peptide_code[tr], y_leak[tr])
        p = clf.predict_proba(peptide_code[te])[:, 1]
        naive_aucs.append(roc_auc_score(y_leak[te], p))
    naive_auc = float(np.mean(naive_aucs))

    # CV-1-style grouped split (groups = peptide cluster)
    gkf = GroupKFold(n_splits=5)
    grouped_aucs = []
    for tr, te in gkf.split(peptide_code, groups=cluster_ids):
        clf = HistGradientBoostingClassifier(categorical_features=[0], **HGB_KWARGS)
        clf.fit(peptide_code[tr], y_leak[tr])
        p = clf.predict_proba(peptide_code[te])[:, 1]
        grouped_aucs.append(roc_auc_score(y_leak[te], p))
    grouped_auc = float(np.mean(grouped_aucs))

    assert naive_auc > 0.9, f"self-test FAILED: naive random split AUROC {naive_auc:.3f} should be ~1 (leaky)"
    assert 0.35 < grouped_auc < 0.65, (
        f"self-test FAILED: peptide-grouped split AUROC {grouped_auc:.3f} should be ~0.5 "
        "(no signal available once peptide identity cannot leak train->test)"
    )

    print(f"Self-test PASSED: recovery rho={rho:.3f} (expected~{abs(expected_rho):.3f}); "
          f"leakage test naive-split AUROC={naive_auc:.3f} (leaky, expected ~1), "
          f"peptide-grouped AUROC={grouped_auc:.3f} (expected ~0.5).")


# =============================================================================
# Figure
# =============================================================================
def plot_summary_figure(all_summaries, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from utils import set_plot_style
        set_plot_style()
    except ImportError:
        import matplotlib as mpl
        mpl.rcParams.update({
            "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
            "font.family": "sans-serif", "font.size": 10,
            "axes.titlesize": 11, "axes.labelsize": 10,
            "axes.spines.top": False, "axes.spines.right": False,
            "xtick.labelsize": 9, "ytick.labelsize": 9,
            "legend.fontsize": 9, "legend.frameon": False, "lines.linewidth": 1.5,
        })

    MM = 1 / 25.4
    COL2_WIDTH = 183 * MM

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, COL2_WIDTH / 2.2), sharey=True)
    for ax, assay in zip(axes, ["IC50", "Kd"]):
        schemes = ["cv1", "cv2"]
        x = np.arange(len(schemes))
        width = 0.25
        for i, model in enumerate(["B1", "B2", "B3"]):
            vals, lo, hi = [], [], []
            for scheme in schemes:
                s = all_summaries.get((assay, scheme))
                if s is None:
                    vals.append(np.nan); lo.append(0); hi.append(0)
                    continue
                m = s["pooled"][model]
                rho = m.get("spearman")
                vals.append(rho if rho is not None else np.nan)
                lo.append((rho - m["spearman_ci_lo"]) if (rho is not None and m["spearman_ci_lo"] is not None) else 0)
                hi.append((m["spearman_ci_hi"] - rho) if (rho is not None and m["spearman_ci_hi"] is not None) else 0)
            color = {"B1": "#CC3311", "B2": "#4477AA", "B3": "#4477AA"}[model]
            alpha = {"B1": 0.9, "B2": 0.45, "B3": 0.9}[model]
            ax.bar(x + (i - 1) * width, vals, width, yerr=[lo, hi], capsize=2,
                   color=color, alpha=alpha, label=model, edgecolor="white", linewidth=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(["peptide-\ngrouped\n(CV-1)", "leave-one-\nallele-out\n(CV-2)"], fontsize=8)
        ax.axhline(0, color="#888888", linewidth=0.6)
        ax.set_title(assay, fontsize=10)
    axes[0].set_ylabel("pooled Spearman ρ\n(log10 affinity)")
    axes[1].legend(["B1: single score", "B2: sequence-only", "B3: sequence+structure"],
                    loc="upper right", fontsize=6.5, frameon=False)
    fig.suptitle("Minimal baseline: pooled Spearman ρ under leakage-controlled CV", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_dir / "affinity_baseline_summary.pdf")
    fig.savefig(out_dir / "affinity_baseline_summary.png", dpi=300)
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=str, default=None,
                     help="Directory to write outputs. Required unless --self-test. "
                          "Never write outputs into the scripts directory.")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--skip-cv3", action="store_true",
                     help="Skip the optional strictest scheme (spec 4.3, CV-3) to save time.")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if args.out_dir is None:
        raise SystemExit("--out-dir is required (never write outputs into the scripts directory)")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    df = load_and_prepare()
    print(f"Loaded {len(df):,} unflagged rows from {METADATA_FN}")

    all_peptides = df[COL_PEPTIDE].unique().tolist()
    pep_to_cluster, cluster_diag = cluster_peptides(all_peptides)
    cluster_diag = cluster_diagnostics_to_dict(cluster_diag, df)
    print("Peptide clustering:", json.dumps(cluster_diag, indent=2))
    df["cluster_id"] = df[COL_PEPTIDE].map(pep_to_cluster)

    pd.DataFrame({"peptide": list(pep_to_cluster.keys()), "cluster_id": list(pep_to_cluster.values())}) \
        .to_csv(out_dir / "peptide_clusters.csv", index=False)

    all_metrics = {"cluster_diagnostics": cluster_diag, "assays": {}}
    fold_rows_all = []
    per_allele_rows_all = []
    all_summaries_for_fig = {}
    kd_sensitivity = None

    for assay_label, meas_type, ceilings, floor in ASSAYS:
        sub = filter_assay(df, meas_type, ceilings, floor)
        print(f"\n=== {assay_label}: {len(sub):,} quantitative, non-censored rows ===")
        assay_out = {"n_rows": int(len(sub)), "schemes": {}}

        schemes = [("cv1", True), ("cv2", False)]
        if not args.skip_cv3:
            schemes.append(("cv3", False))

        for scheme_name, include_allele in schemes:
            reg_oof, n_folds = run_scheme(sub, scheme_name, include_allele, task="regression")
            clf_oof, _ = run_scheme(sub, scheme_name, include_allele, task="classification")
            if len(reg_oof) == 0:
                print(f"  {scheme_name}: no folds (insufficient data), skipped.")
                continue

            reg_summary, reg_per_allele = summarize_scheme(reg_oof, "regression", assay_label, scheme_name)
            clf_summary, clf_per_allele = summarize_scheme(clf_oof, "classification", assay_label, scheme_name)
            reg_summary["n_folds"] = n_folds
            clf_summary["n_folds"] = n_folds

            assay_out["schemes"][scheme_name] = {"regression": reg_summary, "classification": clf_summary}
            all_summaries_for_fig[(assay_label, scheme_name)] = reg_summary

            for m, padf in reg_per_allele.items():
                if len(padf):
                    padf = padf.copy()
                    padf["assay"], padf["scheme"], padf["model"], padf["task"] = assay_label, scheme_name, m, "regression"
                    per_allele_rows_all.append(padf)
            for m, padf in clf_per_allele.items():
                if len(padf):
                    padf = padf.copy()
                    padf["assay"], padf["scheme"], padf["model"], padf["task"] = assay_label, scheme_name, m, "classification"
                    per_allele_rows_all.append(padf)

            reg_oof2 = reg_oof.copy()
            reg_oof2["assay"], reg_oof2["scheme"], reg_oof2["task"] = assay_label, scheme_name, "regression"
            fold_rows_all.append(reg_oof2)
            clf_oof2 = clf_oof.copy()
            clf_oof2["assay"], clf_oof2["scheme"], clf_oof2["task"] = assay_label, scheme_name, "classification"
            fold_rows_all.append(clf_oof2)

            print(f"  {scheme_name} ({n_folds} folds): "
                  f"B3 pooled rho={reg_summary['pooled']['B3']['spearman']:.3f} "
                  f"[{reg_summary['pooled']['B3']['spearman_ci_lo']}, {reg_summary['pooled']['B3']['spearman_ci_hi']}], "
                  f"B2={reg_summary['pooled']['B2']['spearman']:.3f}, B1={reg_summary['pooled']['B1']['spearman']:.3f}, "
                  f"delta(B3-B2)={reg_summary['paired_deltas']['B3_minus_B2']['delta']}")

            if assay_label == "Kd" and scheme_name == "cv1":
                kd_sensitivity = kd_sublabel_sensitivity(reg_oof[reg_oof.columns].copy())

        all_metrics["assays"][assay_label] = assay_out

    if kd_sensitivity is not None:
        all_metrics["kd_sublabel_sensitivity_cv1_B3"] = kd_sensitivity

    with open(out_dir / "affinity_baseline_metrics.json", "w") as f:
        json.dump(all_metrics, f, indent=2)

    if per_allele_rows_all:
        pd.concat(per_allele_rows_all, ignore_index=True).to_csv(out_dir / "affinity_baseline_per_allele.csv", index=False)
    if fold_rows_all:
        fold_cols = ["assay", "scheme", "task", "fold", "held_out_allele", "allele", "cluster_id", "peptide",
                     "y_reg", "y_clf", "pred_B0", "pred_B1", "pred_B2", "pred_B3"]
        combined = pd.concat(fold_rows_all, ignore_index=True)
        combined[fold_cols].to_csv(out_dir / "affinity_baseline_fold_predictions.csv", index=False)
        combined[["assay", "scheme", "task", "fold", "held_out_allele"]].drop_duplicates() \
            .to_csv(out_dir / "fold_assignments.csv", index=False)

    plot_summary_figure(all_summaries_for_fig, out_dir)

    print(f"\nOutputs written to {out_dir}")


if __name__ == "__main__":
    main()
