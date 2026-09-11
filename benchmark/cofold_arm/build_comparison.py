#!/usr/bin/env python
"""Build the combined all_arms_summary.csv and print the cross-arm comparison
report (coverage, paired Wilcoxon vs Rosetta, within-target spread).

Read-only w.r.t. rosetta_arm/results and the per-arm results already produced.
Writes only cofold_arm/results/all_arms_summary.csv.
"""
import csv
import statistics
from pathlib import Path

from scipy.stats import wilcoxon

BASE = Path("/data/p_csb_meiler/huntek1/benchmark")
ARMS = ["af2_reuse", "boltz2", "boltz1", "chai1"]

def load_summary(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {r["pdb_id"]: r for r in rows}

def load_decoys(path):
    from collections import defaultdict
    d = defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                rmsd = float(r["peptide_backbone_rmsd"])
            except (ValueError, KeyError):
                continue
            d[r["pdb_id"]].append(rmsd)
    return d

rosetta = load_summary(BASE / "rosetta_arm/results/arm_summary.csv")
arm_summaries = {arm: load_summary(BASE / f"cofold_arm/{arm}/results/{arm}_summary.csv") for arm in ARMS}
arm_decoys = {arm: load_decoys(BASE / f"cofold_arm/{arm}/results/{arm}_decoys.csv") for arm in ARMS}

all_pdb_ids = sorted(set().union(rosetta.keys(), *[s.keys() for s in arm_summaries.values()]))

# --- write combined table ---
out_path = BASE / "cofold_arm/results/all_arms_summary.csv"
out_path.parent.mkdir(parents=True, exist_ok=True)
fieldnames = ["pdb_id", "peptide_seq"]
for arm in ["rosetta"] + ARMS:
    fieldnames += [f"{arm}_n_scored", f"{arm}_selected", f"{arm}_oracle",
                   f"{arm}_median", f"{arm}_rank_of_selected", f"{arm}_percentile_of_selected"]

with open(out_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames)
    w.writeheader()
    for pdb_id in all_pdb_ids:
        row = {"pdb_id": pdb_id}
        pep = None
        for arm_name, table in [("rosetta", rosetta)] + [(a, arm_summaries[a]) for a in ARMS]:
            r = table.get(pdb_id)
            if r is None:
                for k in ["n_scored", "selected", "oracle", "median", "rank_of_selected", "percentile_of_selected"]:
                    row[f"{arm_name}_{k}"] = ""
                continue
            pep = pep or r["peptide_seq"]
            for k in ["n_scored", "selected", "oracle", "median", "rank_of_selected", "percentile_of_selected"]:
                row[f"{arm_name}_{k}"] = r[k]
        row["peptide_seq"] = pep
        w.writerow(row)

print(f"Wrote {out_path} ({len(all_pdb_ids)} rows)")

# --- coverage ---
print("\n=== COVERAGE ===")
print(f"rosetta: {len(rosetta)} targets")
for arm in ARMS:
    print(f"{arm}: {len(arm_summaries[arm])} targets")

common = set(rosetta.keys())
for arm in ARMS:
    common &= set(arm_summaries[arm].keys())
common = sorted(common)
print(f"\nCommon target set (rosetta + all 4 cofold arms): N = {len(common)}")

# --- per-arm descriptive stats on full arm coverage ---
print("\n=== PER-ARM DESCRIPTIVE STATS (own full coverage) ===")
def stats_block(name, table):
    sel = [float(r["selected"]) for r in table.values() if r["selected"] not in ("", "nan")]
    orc = [float(r["oracle"]) for r in table.values() if r["oracle"] not in ("", "nan")]
    med = [float(r["median"]) for r in table.values() if r["median"] not in ("", "nan")]
    n = len(table)
    under2_sel = sum(1 for x in sel if x < 2.0)
    print(f"{name}: N={n}")
    print(f"  selected: mean={statistics.mean(sel):.3f} median={statistics.median(sel):.3f}  %<2A={100*under2_sel/len(sel):.1f}%")
    print(f"  oracle:   mean={statistics.mean(orc):.3f} median={statistics.median(orc):.3f}")
    print(f"  ensmedian:mean={statistics.mean(med):.3f} median={statistics.median(med):.3f}")

stats_block("rosetta", rosetta)
for arm in ARMS:
    stats_block(arm, arm_summaries[arm])

# --- within-target spread ---
print("\n=== WITHIN-TARGET SPREAD (per target: max-min RMSD across its samples) ===")
def spread_block(name, decoys):
    spreads = []
    for pdb_id, vals in decoys.items():
        if len(vals) >= 2:
            spreads.append(max(vals) - min(vals))
    if spreads:
        print(f"{name}: N_targets={len(spreads)} spread min={min(spreads):.3f} median={statistics.median(spreads):.3f} max={max(spreads):.3f} mean={statistics.mean(spreads):.3f}")

for arm in ARMS:
    spread_block(arm, arm_decoys[arm])

# --- paired comparison vs rosetta on common set ---
print(f"\n=== PAIRED COMPARISON vs ROSETTA (common set, N={len(common)}) ===")
for arm in ARMS:
    diffs = []
    wins_arm = 0
    wins_rosetta = 0
    ties = 0
    for pdb_id in common:
        r_sel = rosetta[pdb_id]["selected"]
        a_sel = arm_summaries[arm][pdb_id]["selected"]
        if r_sel in ("", "nan") or a_sel in ("", "nan"):
            continue
        r_sel = float(r_sel)
        a_sel = float(a_sel)
        diff = a_sel - r_sel  # positive => arm worse (higher RMSD) than rosetta
        diffs.append(diff)
        if a_sel < r_sel:
            wins_arm += 1
        elif r_sel < a_sel:
            wins_rosetta += 1
        else:
            ties += 1
    n = len(diffs)
    mean_diff = statistics.mean(diffs)
    median_diff = statistics.median(diffs)
    try:
        stat, p = wilcoxon(diffs)
    except ValueError as e:
        stat, p = float("nan"), float("nan")
    print(f"{arm} vs rosetta: N={n}  mean(arm-rosetta)={mean_diff:+.3f}  median(arm-rosetta)={median_diff:+.3f}  "
          f"{arm}_wins={wins_arm} ({100*wins_arm/n:.1f}%)  rosetta_wins={wins_rosetta} ({100*wins_rosetta/n:.1f}%)  ties={ties}  "
          f"Wilcoxon stat={stat:.1f} p={p:.4g}")
