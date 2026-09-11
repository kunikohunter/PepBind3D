# Validation pair regeneration

Regenerate FlexPepDock structural ensembles for the 52 validation pairs
using the patched `HLA_db.py` with self-template exclusion enabled
(via the pipeline's `--ignore_epitope_match` flag).

## Why this exists

The released structures for these 52 pairs were threaded onto their own
crystal, because the self-exclusion path in `HLA_db.get_peptide_template` did
not run. The patched `HLA_db.py` fixes it, and this directory re-runs the
validation subset so the analysis measures modelling accuracy rather than
refinement of a self-template. Only pairs whose native PDB is in the local
template database are affected — roughly 0.1% of the 49,268 released
structures.

## Workflow

```
build_validation_fastas.py        # generates per-allele FASTAs from rmsd_per_pair.csv
run_regeneration.sh               # iterates over alleles, calls IEDBTestPipeline.py
```

### Step 1 — Generate FASTAs

```
cd /home/huntek1/main_project/scripts/IEDB_validation/regeneration
python build_validation_fastas.py
```

Outputs per-allele FASTA files into
`/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/fastas/`,
plus an `allele_manifest.csv` summarizing the work.

### Step 2 — Run the regeneration

```
bash run_regeneration.sh
```

For each allele in the manifest, runs `IEDBTestPipeline.py` with the
allele's FASTA and the `--ignore_epitope_match` flag. Output structures
land in
`/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/structures/{allele_dir}/`.

Runtime: roughly 2.5 hours wall time on 28 cores for 52 pairs at 25 decoys
each (~3 min per FlexPepDock refinement).

To run in the background and stay logged out of your shell:
```
nohup bash run_regeneration.sh > regeneration.out 2>&1 &
```

### Step 3 — Rerun the validation notebook

After regeneration completes, point `01_structural_validation.ipynb`'s
`PDB_ROOT` variable at
`/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/structures/`
and rerun. The same 52 pairs will be evaluated; the RMSDs will reflect
threading from non-self templates.

## Prerequisites

- The patched `/home/huntek1/main_project/scripts/HLA_db.py` is in place
  (smoke test on the cluster confirmed `omit=['self']` now filters self-matches).
- `rmsd_per_pair.csv` exists at the expected path from a previous run of
  `01_structural_validation.ipynb`. Without it, there's no list of pairs
  to regenerate.
- The Python environment used has the packages required by the pipeline
  (biopython, pandas, BioPython substitution_matrices, etc.).

## Output structure

```
data/IEDB_data_clean/IEDB_validation/regeneration/
├── allele_manifest.csv
├── fastas/
│   ├── A0201_validation.fasta
│   └── ...
├── logs/
│   ├── A0201.log
│   └── ...
└── structures/
    ├── A0201/
    │   └── (per-peptide subdirs with 25 decoys + score.sc, same layout
    │        as the released pdb/ directory)
    └── ...
```
