# Re-dock the 24 newly crystal-matched pairs (ACCRE)

## Why

`analysis/crystal_match.py` finds 76 pairs in the merged release with a matching
crystal structure in the local MHC template database: the original 52 plus **24
new ones from the v2 batch**.

The released ensembles for all of them are self-templated, threaded onto their
own crystal, because `HLA_db.get_peptide_template` ranks same-length peptides
from the whole gene by BLOSUM62 similarity (so an identical peptide wins) and
self-exclusion runs only under `--ignore_epitope_match`, which production did
not pass. Their RMSD to that crystal therefore measures refinement of a
self-template, not modelling accuracy.

The original 52 already have leakage-free ensembles in
`IEDB_validation/regeneration/` (re-docked with the flag). **These 24 do not.**
Re-docking them is what lets Validation 1 report 76 pairs instead of 52.

Do **not** re-dock the original 52, their existing regeneration ensembles were
produced by this same recipe and re-running adds nothing.

## Cost

24 pairs x 25 decoys = 600 decoys, 1 CPU per task, ~2.5 h walltime each as a
task array: **≤ 60 CPU-hours**, no GPU. Output is ~600 PDBs plus scorefiles,
well under 1 GB.

## Rosetta build: use the v2 one

These 24 pairs come from the v2 batch, generated with
**2024.09+release.06b3cf8**. Use that build, not the 3.15 tree that
`run_regeneration.sh` hardcodes for Tungsten. Mixing a third Rosetta version
into a number the paper reports is what we are trying to avoid.

## Inputs to copy to ACCRE

Already generated on Tungsten:

```
.../IEDB_data_clean/IEDB_validation/crystal_match_v2/redock_24_pairs.csv
.../IEDB_data_clean/IEDB_validation/regeneration_v2/fastas/*_validation_v2.fasta   # 10 files
.../IEDB_data_clean/IEDB_validation/regeneration_v2/allele_manifest_validation_v2.csv
```

10 alleles: A*02:01 (12 peptides), A*11:01 (3), A*24:02 (2), and one each for
A*02:03, A*02:06, A*02:07, A*03:01, A*26:01, B*07:02, B*39:01.

```bash
# from Tungsten
rsync -avP \
  $PEPBIND3D_DATA/IEDB_validation/regeneration_v2/ \
  $PEPBIND3D_DATA/IEDB_validation/crystal_match_v2/redock_24_pairs.csv \
  $CLUSTER_LOGIN:$PEPBIND3D_CLUSTER/main_project/data/redock_24/
```

## Steps on ACCRE

Two stages, run separately, same as the existing threading/docking arrays.

**1. Threading, the `--ignore_epitope_match` flag is the whole point.**

Per allele, following `run_regeneration.sh:60-68`:

```bash
python $REPO/../IEDBTestPipeline_ACCRE.py \
    --IEDBquery skip \
    --buildFasta <allele>_validation_v2.fasta 0 \
    --setAllele "<allele in A*02:01 form>" \
    --rosetta <ROSETTA_2024.09_MAIN> \
    --threads <n> \
    --ignore_epitope_match
```

Verify the patched `HLA_db.py` (the one whose `get_peptide_template` honours
`omit=["self"]`) is the copy on ACCRE's `PYTHONPATH`. Smoke-test on one allele
first and confirm from the log that the chosen template is **not** the pair's
own crystal, grep the `Running Rosetta with template` line, or the
`{peptide}_ROSETTA.log` `-s .../templates/{PDB}.pdb` argument, against
`redock_24_pairs.csv`'s `matched_pdb_id` column. If any pair still threads onto
its own `matched_pdb_id`, stop: the flag is not taking effect and the re-dock is
pointless.

**2. FlexPepDock refinement, 25 decoys, production settings.**

```
-pep_refine -nstruct 25 -ex1 -ex2aro
```

(with the prepack step first, as in the production pipeline). Then check every
pair produced exactly 25 decoys and a `score.sc` whose `description` column
matches the PDB stems, the release conversion hard-fails otherwise.

## Return to Tungsten

```bash
rsync -avP \
  $CLUSTER_LOGIN:$PEPBIND3D_CLUSTER/main_project/data/redock_24/output/ \
  $PEPBIND3D_DATA/IEDB_validation/regeneration_v2/pdb/
```

Expected layout, matching the v1 regeneration tree so the analysis can read
both: `regeneration_v2/pdb/{allele_dir}/{peptide}/` containing
`{peptide}_input_{0001..0025}.pdb` and `score.sc`.

## Then, on Tungsten

Score the re-docked ensembles and fold them into Validation 1. These are
**PDB-tree** inputs (like the v1 regeneration), not silent files, so they go
through the notebook-01 code path rather than `analysis/crystal_rmsd.py`, which
reads release silents and reports upper bounds only.

Sanity check on arrival: the re-docked RMSDs should be **worse** than the
released-structure upper bounds already measured for these same 24 pairs
(median 0.96 Å best-score, in
`IEDB_validation/crystal_rmsd_v2/crystal_rmsd_per_pair.csv`). On the original 52
removing the leakage cost ~0.24 Å, and on that run ~0.33 Å, so expect
roughly 1.2-1.3 Å. If the re-docked numbers come back *better*, the flag did not
work, do not report them.
