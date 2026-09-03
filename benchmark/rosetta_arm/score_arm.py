"""Score docked pMHC decoys for the full arm run and record ensemble rank.

Extends score_pilot.py (kept intact) with, per target:
  rank_of_selected        1-indexed rank (ascending peptide_backbone_rmsd) of
                           the decoy Rosetta itself ranks best by total_score.
                           1 = Rosetta picked the most accurate decoy in the
                           ensemble; n_scored = it picked the worst.
  percentile_of_selected   (rank_of_selected - 1) / (n_scored - 1), in [0, 1],
                           so targets with fewer than 25 scored decoys stay
                           comparable to those with a full 25. NaN if
                           n_scored == 1.

Also keeps selected / oracle / median / n_scored per target, and writes a
tidy long-format CSV with all per-decoy RMSDs and total_scores so ensemble
distributions can be analyzed later without re-scoring.

Outputs (written to /data/p_csb_meiler/huntek1/benchmark/rosetta_arm/results/):
  arm_summary.csv   one row per target: selected/oracle/median/rank/percentile/n
  arm_decoys.csv    one row per (target, decoy): rmsd, total_score
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, "/data/p_csb_meiler/huntek1/benchmark/metric")
from pmhc_rmsd import score_prediction

ARM = Path("/data/p_csb_meiler/huntek1/benchmark/rosetta_arm/arm")
REFS = Path("/data/p_csb_meiler/huntek1/benchmark/refs/cif")
MANIFEST = Path("/data/p_csb_meiler/huntek1/benchmark/refs/reference_manifest.csv")
RESULTS = Path("/data/p_csb_meiler/huntek1/benchmark/rosetta_arm/results")
RESULTS.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = RESULTS / "arm_summary.csv"
DECOYS_CSV = RESULTS / "arm_decoys.csv"

peptides = {r["pdb_id"]: r["peptide_seq"] for r in csv.DictReader(open(MANIFEST))}


def read_scores(score_sc):
    """Return {decoy_tag: total_score} from a Rosetta score.sc."""
    out, header = {}, None
    for line in open(score_sc):
        if not line.startswith("SCORE:"):
            continue
        parts = line.split()[1:]
        if header is None:
            header = parts
            continue
        if len(parts) != len(header):
            continue
        row = dict(zip(header, parts))
        try:
            out[row["description"]] = float(row["total_score"])
        except (KeyError, ValueError):
            pass
    return out


summary_rows = []
decoy_rows = []

print(f"{'target':8s} {'pep':11s} {'n':>3s} {'selected':>9s} {'oracle':>8s} {'median':>8s} {'rank':>5s} {'pct':>6s}")
print("-" * 72)

for tdir in sorted(ARM.iterdir()):
    if not tdir.is_dir():
        continue
    tid = tdir.name
    pep = peptides.get(tid)
    ref = REFS / f"{tid}.cif"
    if not pep or not ref.exists():
        print(f"{tid:8s} SKIP (no peptide or reference)")
        continue

    decoys = sorted(tdir.glob("docking/*_[0-9][0-9][0-9][0-9].pdb"))
    if not decoys:
        print(f"{tid:8s} SKIP (no decoys)")
        continue

    score_sc = tdir / "docking" / "score.sc"
    if not score_sc.exists():
        print(f"{tid:8s} SKIP (no score.sc)")
        continue
    scores = read_scores(score_sc)

    rms = {}
    for d in decoys:
        try:
            rms[d.stem] = score_prediction(d, ref, pep, copy_policy="best")["peptide_backbone_rmsd"]
        except Exception as e:
            print(f"   ! {tid} {d.name}: {type(e).__name__}: {e}", file=sys.stderr)

    if not rms:
        print(f"{tid:8s} SKIP (all decoys failed to score)")
        continue

    for tag, rmsd in rms.items():
        decoy_rows.append({
            "pdb_id": tid,
            "peptide_seq": pep,
            "decoy": tag,
            "peptide_backbone_rmsd": rmsd,
            "total_score": scores.get(tag, ""),
        })

    # Rank decoys by ascending RMSD (1 = most accurate) among those we scored.
    rmsd_ranked = sorted(rms, key=lambda t: rms[t])
    rmsd_rank = {tag: i + 1 for i, tag in enumerate(rmsd_ranked)}
    n_scored = len(rms)

    # Rosetta's own selection: lowest total_score among decoys we could score.
    score_ranked = [t for t in sorted(scores, key=scores.get) if t in rms]
    if score_ranked:
        selected_tag = score_ranked[0]
        sel = rms[selected_tag]
        rank_of_selected = rmsd_rank[selected_tag]
        percentile_of_selected = (
            (rank_of_selected - 1) / (n_scored - 1) if n_scored > 1 else float("nan")
        )
    else:
        sel = float("nan")
        rank_of_selected = None
        percentile_of_selected = float("nan")
        print(f"         (warning: no score.sc/RMSD tag overlap; 'selected' unavailable)")

    vals = sorted(rms.values())
    oracle = min(vals)
    med = vals[len(vals) // 2]

    summary_rows.append({
        "pdb_id": tid,
        "peptide_seq": pep,
        "n_scored": n_scored,
        "selected": sel,
        "oracle": oracle,
        "median": med,
        "rank_of_selected": rank_of_selected if rank_of_selected is not None else "",
        "percentile_of_selected": percentile_of_selected,
    })

    rank_str = str(rank_of_selected) if rank_of_selected is not None else "NA"
    pct_str = f"{percentile_of_selected:.3f}" if percentile_of_selected == percentile_of_selected else "NA"
    print(f"{tid:8s} {pep:11s} {n_scored:3d} {sel:9.3f} {oracle:8.3f} {med:8.3f} {rank_str:>5s} {pct_str:>6s}")

with open(SUMMARY_CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=[
        "pdb_id", "peptide_seq", "n_scored", "selected", "oracle", "median",
        "rank_of_selected", "percentile_of_selected",
    ])
    w.writeheader()
    w.writerows(summary_rows)

with open(DECOYS_CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=[
        "pdb_id", "peptide_seq", "decoy", "peptide_backbone_rmsd", "total_score",
    ])
    w.writeheader()
    w.writerows(decoy_rows)

print(f"\nWrote {len(summary_rows)} target rows to {SUMMARY_CSV}")
print(f"Wrote {len(decoy_rows)} decoy rows to {DECOYS_CSV}")
