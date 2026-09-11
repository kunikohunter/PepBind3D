"""Score docked pMHC decoys against their experimental reference structures.

Reports three numbers per target:
  selected  - RMSD of the decoy Rosetta itself ranks best (lowest total_score).
              This is the only prospectively honest number: it is what the
              protocol returns when the answer is unknown.
  oracle    - RMSD of the best decoy in the ensemble. Unreachable in practice;
              it bounds what the ensemble contains.
  median    - ensemble centre, for spread.
The selected-vs-oracle gap is the cost of scoring error, and the oracle-vs-
median gap is what generating an ensemble (rather than one structure) buys.
"""
import csv, sys, re
from pathlib import Path

sys.path.insert(0, "/data/p_csb_meiler/huntek1/benchmark/metric")
from pmhc_rmsd import score_prediction

PILOT = Path("/data/p_csb_meiler/huntek1/benchmark/rosetta_arm/pilot")
REFS = Path("/data/p_csb_meiler/huntek1/benchmark/refs/cif")
MANIFEST = Path("/data/p_csb_meiler/huntek1/benchmark/refs/reference_manifest.csv")

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


print(f"{'target':8s} {'pep':11s} {'n':>3s} {'selected':>9s} {'oracle':>8s} {'median':>8s}")
print("-" * 56)
for tdir in sorted(PILOT.iterdir()):
    if not tdir.is_dir():
        continue
    tid = tdir.name
    pep = peptides.get(tid)
    ref = REFS / f"{tid}.cif"
    if not pep or not ref.exists():
        print(f"{tid:8s} SKIP (no peptide or reference)")
        continue

    decoys = sorted(tdir.glob("docking/*_0[0-9][0-9][0-9].pdb"))
    if not decoys:
        print(f"{tid:8s} SKIP (no decoys)")
        continue

    scores = read_scores(tdir / "docking" / "score.sc")
    rms = {}
    for d in decoys:
        try:
            res = score_prediction(d, ref, pep, copy_policy="best")
        except Exception as e:
            print(f"   ! {tid} {d.name}: {type(e).__name__}: {e}", file=sys.stderr)
            continue
        rmsd = res["peptide_backbone_rmsd"]
        assert rmsd <= 10.0, (
            f"{tid} {d.name}: peptide_backbone_rmsd={rmsd:.3f} exceeds 10 A "
            f"sanity ceiling -- treat as a correspondence failure, not a bad model"
        )
        rms[d.stem] = rmsd
        for w in res.get("warnings", []) or []:
            print(f"   ! {tid} {d.name}: {w}", file=sys.stderr)

    if not rms:
        print(f"{tid:8s} SKIP (all decoys failed to score)")
        continue

    ranked = [t for t in sorted(scores, key=scores.get) if t in rms]
    sel = rms[ranked[0]] if ranked else float("nan")
    vals = sorted(rms.values())
    med = vals[len(vals) // 2]
    print(f"{tid:8s} {pep:11s} {len(rms):3d} {sel:9.3f} {min(vals):8.3f} {med:8.3f}")
    if not ranked:
        print(f"         (warning: no score.sc/RMSD tag overlap; 'selected' unavailable)")
