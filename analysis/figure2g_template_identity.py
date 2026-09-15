"""
Figure 2G, rebuilt from the threading logs: best-score RMSD vs the identity of
the template ACTUALLY used.

The submitted Figure 2G reports Spearman rho = -0.42 (p = 0.002, n = 50)
between best-score RMSD and threading-template identity. Those identities came
from `best_template_for_pair` in 01_structural_validation.ipynb (cell 19), a
reimplementation of the selector that filters candidates to the SAME ALLELE:

    candidates = mhc_db_info[(mhc_db_info['MHC_Allele'] == allele) & ...]

The real selector (HLA_db.MHCdatabase.get_peptide_template) pools every
same-length peptide from every allele of the same GENE and ranks by BLOSUM62
similarity. So the notebook drew from a smaller candidate pool and could name a
template that was never used: it disagrees with the run logs on 9 of the 50
pairs, and the correlation is built on those values.

This script takes the template from the source of truth instead -- the
`-s .../templates/{PDB}.pdb` argument recorded in each pair's
`{peptide}_ROSETTA.log` from the regeneration run -- and recomputes the
correlation. Identity is computed exactly as the notebook did (positional
matches / peptide length * 100), so the only thing that changes is WHICH
template each identity describes.

RMSDs come from the regeneration analysis (01_structural_regen/rmsd_per_pair.csv),
which is the leakage-free measurement: the regeneration run passed
--ignore_epitope_match, and these logs confirm 0 of 52 pairs used a self
template.

Usage:
    python3 figure2g_template_identity.py --out-dir <dir>
    python3 figure2g_template_identity.py --self-test
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT, MHC_DB_ROOT  # noqa: E402

BASE = DATA_ROOT
LOG_ROOT = BASE / "IEDB_validation" / "regeneration" / "fastas" / "output"
RMSD_CSV = BASE / "IEDB_validation" / "01_structural_regen" / "rmsd_per_pair.csv"
MHC_DB = MHC_DB_ROOT / "database.info"

# "-s <MHC_DB_ROOT>/templates/4NQX.pdb"
TEMPLATE_RE = re.compile(r"templates/([0-9A-Za-z]{4})\.pdb")
# "A0201_validation_batch1" -> "A0201"
DIRNAME_RE = re.compile(r"^([A-Z][0-9]{4})_")


def allele_to_dir(allele: str) -> str:
    s = allele[4:] if allele.startswith("HLA-") else allele
    return s.replace("*", "").replace(":", "")


def percent_identity(query: str, template: str) -> float:
    """Positional identity, as the notebook computed it. The selector requires
    equal length, so no alignment is needed; unequal lengths mean something is
    wrong upstream and the caller is told."""
    if len(query) != len(template):
        raise ValueError(f"length mismatch: {query} ({len(query)}) vs "
                         f"{template} ({len(template)})")
    return sum(a == b for a, b in zip(query, template)) / len(query) * 100.0


def templates_from_logs(log_root=LOG_ROOT):
    """(allele_dir, peptide) -> template PDB ID actually passed to Rosetta.

    Keyed on allele AND peptide: two validation pairs can share a peptide on
    different alleles, so keying on peptide alone would silently mis-assign."""
    rows = []
    for log in sorted(log_root.rglob("*_ROSETTA.log")):
        m = TEMPLATE_RE.search(log.read_text(errors="ignore"))
        if not m:
            continue
        d = DIRNAME_RE.match(log.parent.parent.name)
        if not d:
            continue
        rows.append({"allele_dir": d.group(1), "peptide": log.parent.name,
                     "template_used": m.group(1).upper()})
    df = pd.DataFrame(rows)
    dup = df.duplicated(subset=["allele_dir", "peptide"], keep=False)
    if dup.any():
        raise SystemExit("two logs for the same (allele, peptide):\n"
                         f"{df[dup].to_string(index=False)}")
    return df


# The 24 newly matched pairs were re-docked separately and their logs arrived
# in a different layout: <allele_compact>/<peptide>_ROSETTA.log rather than
# <allele>_batch/<peptide>/<peptide>_ROSETTA.log. Namespaced by allele
# deliberately -- FLPSDFFPSV appears on three A*02 subtypes with three different
# crystals, so flat filenames would have collided.
LOG_ROOT_24 = BASE / "results_incoming" / "rosetta_logs_24"
RMSD_CSV_24 = BASE / "IEDB_validation" / "crystal_rmsd_24" / "crystal_rmsd_per_pair.csv"


def templates_from_flat_logs(log_root=LOG_ROOT_24):
    """(allele_dir, peptide) -> template PDB, for <allele>/<peptide>_ROSETTA.log."""
    rows = []
    for log in sorted(log_root.rglob("*_ROSETTA.log")):
        m = TEMPLATE_RE.search(log.read_text(errors="ignore"))
        if not m:
            continue
        rows.append({"allele_dir": log.parent.name,
                     "peptide": log.stem.replace("_ROSETTA", ""),
                     "template_used": m.group(1).upper()})
    df = pd.DataFrame(rows)
    dup = df.duplicated(subset=["allele_dir", "peptide"], keep=False)
    if dup.any():
        raise SystemExit(f"two logs for the same (allele, peptide):\n"
                         f"{df[dup].to_string(index=False)}")
    return df


def build(rmsd_csv=RMSD_CSV, log_root=LOG_ROOT, mhc_db=MHC_DB):
    rmsd = pd.read_csv(rmsd_csv)
    rmsd["allele_dir"] = rmsd["allele"].map(allele_to_dir)
    # Both log trees. rmsd_per_pair.csv carries all 76 crystal-matched pairs
    # (notebook 01 reads both re-docked trees), so reading only the first tree
    # leaves the 24 without a template and silently reports the old 52-pair rho.
    logs = pd.concat([templates_from_logs(log_root),
                      templates_from_flat_logs(LOG_ROOT_24)], ignore_index=True)
    db = pd.read_csv(mhc_db)
    pdb_to_pep = (db.assign(PDB_ID=db["PDB_ID"].astype(str).str.upper())
                    .drop_duplicates(subset=["PDB_ID"])
                    .set_index("PDB_ID")["Epitope_Description"].astype(str).to_dict())

    d = rmsd.merge(logs, on=["allele_dir", "peptide"], how="left")
    d["template_peptide_used"] = d["template_used"].map(pdb_to_pep)

    idents, problems = [], []
    for r in d.itertuples(index=False):
        tp = getattr(r, "template_peptide_used")
        if not isinstance(tp, str) or not isinstance(r.template_used, str):
            idents.append(np.nan)
            problems.append({"allele": r.allele, "peptide": r.peptide,
                             "template_used": r.template_used,
                             "reason": "template not recovered or not in database.info"})
            continue
        try:
            idents.append(percent_identity(r.peptide, tp))
        except ValueError as e:
            idents.append(np.nan)
            problems.append({"allele": r.allele, "peptide": r.peptide,
                             "template_used": r.template_used, "reason": str(e)})
    d["template_identity_used"] = idents

    # is the template its own crystal? (must be none -- the run used the flag)
    d["is_self_template"] = (d["template_used"].astype(str).str.upper()
                             == d["matched_pdb_id"].astype(str).str.upper())
    return d, pd.DataFrame(problems)


def correlate(d, rmsd_col="rmsd_best_score", ident_col="template_identity_used"):
    s = d.dropna(subset=[rmsd_col, ident_col])
    rho = spearmanr(s[ident_col], s[rmsd_col])
    r = pearsonr(s[ident_col], s[rmsd_col])
    return {"n": int(len(s)), "spearman_rho": float(rho.statistic),
            "spearman_p": float(rho.pvalue),
            "pearson_r": float(r.statistic), "pearson_p": float(r.pvalue),
            "identity_min": float(s[ident_col].min()),
            "identity_median": float(s[ident_col].median()),
            "identity_max": float(s[ident_col].max())}


def self_test():
    """Known answers for the two things this script computes."""
    # identity: exact, one mismatch, all different
    assert percent_identity("ABCDE", "ABCDE") == 100.0
    assert percent_identity("ABCDE", "ABCDX") == 80.0
    assert percent_identity("AAAAA", "CCCCC") == 0.0
    try:
        percent_identity("ABC", "ABCD")
        raise AssertionError("unequal lengths must raise")
    except ValueError:
        pass

    # a perfectly monotone decreasing relation must give rho = -1 exactly
    d = pd.DataFrame({"template_identity_used": [10., 20., 30., 40., 50.],
                      "rmsd_best_score": [5., 4., 3., 2., 1.]})
    res = correlate(d)
    assert abs(res["spearman_rho"] + 1.0) < 1e-12, res
    assert res["n"] == 5

    # NaNs must be dropped, not propagated
    d2 = pd.concat([d, pd.DataFrame({"template_identity_used": [np.nan],
                                     "rmsd_best_score": [9.]})], ignore_index=True)
    assert correlate(d2)["n"] == 5

    # log parsing keys on (allele, peptide), not peptide alone
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        for allele, pep, pdb in [("A0201", "SAMEPEP", "1AAA"),
                                 ("B0702", "SAMEPEP", "2BBB")]:
            p = root / f"{allele}_validation_batch1" / pep
            p.mkdir(parents=True)
            (p / f"{pep}_ROSETTA.log").write_text(
                f"blah -s /x/templates/{pdb}.pdb -out:whatever\n")
        got = templates_from_logs(root).set_index("allele_dir")["template_used"].to_dict()
        assert got == {"A0201": "1AAA", "B0702": "2BBB"}, got
    print("Self-test PASSED: identity is exact and rejects unequal lengths, a "
          "monotone relation gives rho = -1, NaNs are dropped, and one peptide on "
          "two alleles resolves to two different templates.")


def plot(d, out, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mpl = matplotlib
    mpl.rcParams.update({"savefig.dpi": 300, "savefig.bbox": "tight",
                         "font.family": "sans-serif", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    MM = 1 / 25.4
    s = d.dropna(subset=["rmsd_best_score", "template_identity_used"])
    fig, ax = plt.subplots(figsize=(89 * MM, 75 * MM))
    ax.scatter(s["template_identity_used"], s["rmsd_best_score"],
               facecolors="none", edgecolors="#4477AA", linewidth=1.1, s=34)
    m, b = np.polyfit(s["template_identity_used"], s["rmsd_best_score"], 1)
    xs = np.linspace(s["template_identity_used"].min(), s["template_identity_used"].max(), 50)
    ax.plot(xs, m * xs + b, color="#333333", linewidth=1.3)
    ax.axhline(2.0, color="#CC3311", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_xlabel("threading-template identity to native peptide (%)")
    ax.set_ylabel("best-score peptide backbone RMSD (Å)")
    ax.set_title(f"Spearman ρ = {res['spearman_rho']:.2f} "
                 f"(p = {res['spearman_p']:.3g}, n = {res['n']})", fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure2g_template_identity.{ext}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test(); return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    self_test(); print()

    d, problems = build()
    print(f"{len(d)} validation pairs; template recovered from logs for "
          f"{int(d['template_used'].notna().sum())}")
    n_self = int(d["is_self_template"].sum())
    print(f"pairs threaded onto their own crystal: {n_self} "
          f"(must be 0 -- the regeneration run passed --ignore_epitope_match)")
    if n_self:
        raise SystemExit("a regeneration pair used a self template; the RMSDs in "
                         "rmsd_per_pair.csv are not leakage-free after all")
    if len(problems):
        print(f"\n{len(problems)} pair(s) with no usable identity:")
        print(problems.to_string(index=False))

    # how far off was the notebook?
    if "template_pdb" in d.columns:
        both = d.dropna(subset=["template_used", "template_pdb"])
        agree = int((both["template_used"].str.upper()
                     == both["template_pdb"].astype(str).str.upper()).sum())
        print(f"\nnotebook's template vs the log's: agree {agree}/{len(both)}")
        diff = both[both["template_used"].str.upper()
                    != both["template_pdb"].astype(str).str.upper()]
        if len(diff):
            print("disagreements (notebook filtered to the same allele; the real "
                  "selector pools the whole gene):")
            print(diff[["allele", "peptide", "template_pdb", "template_identity",
                        "template_used", "template_identity_used"]]
                  .to_string(index=False))

    res = correlate(d)
    res_nb = (correlate(d, ident_col="template_identity")
              if "template_identity" in d.columns else None)
    print(f"\n=== Figure 2G, templates from the logs ===")
    print(f"  Spearman rho = {res['spearman_rho']:.3f} (p = {res['spearman_p']:.4g}, "
          f"n = {res['n']})")
    print(f"  identity range {res['identity_min']:.1f}-{res['identity_max']:.1f}%, "
          f"median {res['identity_median']:.1f}%")
    if res_nb:
        print(f"  for comparison, notebook identities: rho = {res_nb['spearman_rho']:.3f} "
              f"(p = {res_nb['spearman_p']:.4g}, n = {res_nb['n']}) "
              f"-- the submitted figure reports -0.42, p = 0.002, n = 50")

    d.to_csv(out / "figure2g_template_identity.csv", index=False)
    if len(problems):
        problems.to_csv(out / "figure2g_problems.csv", index=False)
    with open(out / "figure2g_template_identity.json", "w") as f:
        json.dump({"from_logs": res, "from_notebook_reimplementation": res_nb}, f, indent=2)
    plot(d, out, res)
    print(f"\nwrote figure2g_template_identity.{{pdf,png,csv,json}} to {out}")


if __name__ == "__main__":
    main()
