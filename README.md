# IEDB pHLA dataset — validation analyses

Validation analyses for the released pHLA dataset
(https://huggingface.co/datasets/kunikohunter/IEDB_pHLA_binding_data).

## Workflow

Two notebooks, run in order:

1. `01_structural_validation.ipynb` — peptide-backbone RMSD between FlexPepDock
   decoys and experimental crystal structures, for the subset of
   peptide–allele pairs with an `assay_pdb_id` annotation in the metadata.
   Reports three metrics per pair: best-score decoy RMSD, mean RMSD across
   the top-5 by score, and lowest RMSD across all 25 decoys. Generates Figure 2.
2. `02_score_affinity_validation.ipynb` — Spearman correlation of Rosetta
   summary scores (`rosetta_best_score`, `rosetta_mean_score`) against
   experimental log(IC50) and log(Kd). Reports per-allele and pooled
   correlations. Generates Figure 3. Reads only metadata.csv; does not
   touch structure files.

Shared utilities are in the `utils/` package.

## Paths

These are hardcoded at the top of each notebook. Edit if your layout differs.

- Released dataset: `/home/huntek1/main_project/data/IEDB_data_clean/huggingface/`
- Local PDB store (one folder per peptide): `/home/huntek1/main_project/data/IEDB_data_clean/pdb/{allele}/{peptide}/`
- Template PDB database: `/home/huntek1/Data/MHC_database/templates/`
- Validation outputs: `/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/`

## Decoy file convention

Each `pdb/{allele}/{peptide}/` folder contains 25 decoys named
`{peptide}_input_{0001..0025}.pdb` and a `score.sc` Rosetta scorefile.
The `description` column in `score.sc` matches the PDB filename stem.

## Chain conventions

- **Modeled structures**: chain A = MHC binding cleft (α₁/α₂ domains,
  ~180 residues); chain B = peptide.
- **Experimental crystal structures**: chain IDs vary across PDB entries.
  The peptide chain is identified by sequence match (or by length, 7–15
  residues, as fallback). The MHC heavy chain is identified as the
  remaining chain in the 170–290 residue range.

## Reproducibility note for users of the HuggingFace release

The `pdb/` directory is not part of the public release (it would roughly
double the dataset size). To reproduce these analyses from the public
release, first extract the silent files into per-decoy PDBs:

```bash
# Per peptide-allele pair, from inside structures/{allele}/
extract_pdbs.linuxgccrelease -in:file:silent {peptide}.silent
```

The score.sc file is recoverable from the silent file headers if not
included separately.

## Dependencies

Standard scientific Python: `pandas`, `numpy`, `scipy`, `matplotlib`,
`biopython`, `tqdm`. No Rosetta/PyRosetta dependency for analysis (only
required for the upstream `extract_pdbs` step above).
