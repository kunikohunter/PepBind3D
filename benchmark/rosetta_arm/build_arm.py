"""Build the working tree and joblist for the full 174-target Rosetta arm.

Reuses the pilot's exact layout (per target: fasta/ + docking/, one FASTA
containing the reference peptide sequence) and threading invocation
(IEDBTestPipeline_ACCRE.py --IEDBquery skip --setAllele ... --buildFasta
... 0 --rosetta <shim>).

Targets are the 174 rows in resolved_alleles.csv with usable_for_rosetta ==
"yes". A target is SKIPPED (not built, not queued) if:
  - its resolved_allele has no receptor PDB in default_receptor/, or
  - its peptide_seq contains "X" (NCAA placeholder) with no params file
    available (the pilot pipeline never exercised the --params path).

Writes:
  arm/<pdb_id>/fasta/<pdb_id>.fasta
  arm/<pdb_id>/docking/            (created empty; FlexPepDock fills it)
  arm_joblist.txt   one record per queued target, '|'-delimited:
      pdb_id|resolved_allele|peptide_seq|target_dir
  arm_skipped.csv   pdb_id,resolved_allele,reason
"""
import csv
import os
from pathlib import Path

BENCH = Path("/data/p_csb_meiler/huntek1/benchmark")
ALLELES_CSV = BENCH / "alleles/resolved_alleles.csv"
MANIFEST_CSV = BENCH / "refs/reference_manifest.csv"
RECEPTOR_DIR = Path("/home/huntek1/Data/MHC_database/default_receptor")
ARM_ROOT = BENCH / "rosetta_arm/arm"
JOBLIST = BENCH / "rosetta_arm/arm_joblist.txt"
SKIPPED = BENCH / "rosetta_arm/arm_skipped.csv"

alleles = [r for r in csv.DictReader(open(ALLELES_CSV)) if r["usable_for_rosetta"] == "yes"]
manifest = {r["pdb_id"]: r for r in csv.DictReader(open(MANIFEST_CSV))}
have_receptor = set(os.listdir(RECEPTOR_DIR))

queued = []
skipped = []

for row in sorted(alleles, key=lambda r: r["pdb_id"]):
    pdb_id = row["pdb_id"]
    allele = row["resolved_allele"]
    man = manifest.get(pdb_id)
    if man is None:
        skipped.append((pdb_id, allele, "no reference_manifest.csv row"))
        continue
    pep = man["peptide_seq"]
    if f"{allele}.pdb" not in have_receptor:
        reason = "no receptor PDB in default_receptor/"
        if allele == "B*27:13":
            reason += " (receptor rebuilding, job 13779102; do not substitute/approximate)"
        skipped.append((pdb_id, allele, reason))
        continue
    if "X" in pep:
        skipped.append((pdb_id, allele, f"peptide '{pep}' contains NCAA placeholder 'X', no params file"))
        continue

    tdir = ARM_ROOT / pdb_id
    fasta_dir = tdir / "fasta"
    dock_dir = tdir / "docking"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    dock_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = fasta_dir / f"{pdb_id}.fasta"
    fasta_path.write_text(f">{pep}\n{pep}\n")

    queued.append((pdb_id, allele, pep, str(tdir)))

with open(JOBLIST, "w") as fh:
    for pdb_id, allele, pep, tdir in queued:
        fh.write(f"{pdb_id}|{allele}|{pep}|{tdir}\n")

with open(SKIPPED, "w") as fh:
    w = csv.writer(fh)
    w.writerow(["pdb_id", "resolved_allele", "reason"])
    w.writerows(skipped)

print(f"Queued: {len(queued)}")
print(f"Skipped: {len(skipped)}")
for pdb_id, allele, reason in skipped:
    print(f"  SKIP {pdb_id} ({allele}): {reason}")
