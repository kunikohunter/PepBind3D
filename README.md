# PepBind3D — curation and validation code

Code for the data descriptor *PepBind3D: Curated peptide-HLA class I binding
affinities with Rosetta structural ensembles*. This repository contains the
notebooks that produce the manuscript figures and supplementary tables, the
benchmark comparing Rosetta docking against co-folding models, and the
analysis scripts.

The released dataset is on HuggingFace:
https://huggingface.co/datasets/kunikohunter/PepBind3D

The structure-generation pipeline itself is **not** in this repository. It was
published with a prior study and is available there (see below).

## Curation and structure generation

These scripts are published with Bloodworth N, Chen W, Hunter K, Patrick D,
et al. *Posttranslationally modified self-peptides promote hypertension in
mouse models.* J Clin Invest. 2024;134(16):e174374.
doi:10.1172/JCI174374 — code at
https://github.com/meilerlab/discovery-self-peptides-hypertension (`code/`).

This dataset was generated with an adapted copy of that pipeline. Relative to
the published version, the adaptation adds SLURM batch-array execution
(`--batch_index`, `--slurm_setup`, `--threads`, `create_batch_list`,
`safe_thread_all`) and column-standardization helpers, and points at a newer
Rosetta tree; `formats.py` is unchanged. The generation logic — threading
order, template selection, docking protocol — is the same. Structures in this
release stamp `REMARK 220 VERSION 2024.09+release.06b3cf8`.

Note that threading-template selection in the published `HLA_db.py` computes
`omit_self` but does not apply it, so a pair whose own crystal is in the
template database is threaded onto itself. This affects only pairs with a
native PDB in the local database (~0.1% of the release); `regeneration/`
re-runs the 52 validation pairs with self-exclusion enabled so that the
structural validation measures modeling accuracy rather than refinement of a
self-template.

- `IEDBTestPipeline.py` — curation pipeline: reads the IEDB MHC ligand bulk
  download, filters to HLA-A/HLA-B with quantitative IC50/Kd measurements,
  deduplicates per allele, and prepares threaded starting models. Threading
  proceeds SimpleThreadingMover → NCAA substitution → trim → add receptor →
  FastRelax (5 repeats) → FlexPepDock prepack, producing the relaxed starting
  model used as input to docking.
- `HLA_db.py` — builds and queries the local MHC template database used for
  threading-template selection (same-allele, same-length, best BLOSUM62).

## Validation and figure notebooks

Run in order; each reads `metadata.csv` and, where noted, the per-pair Rosetta
score files.

1. `01_structural_validation.ipynb` — peptide-backbone RMSD between FlexPepDock
   decoys and matched experimental crystal structures, for the 52 pairs with a
   PDB match. Reports, per pair, the RMSD of the best decoy (lowest I_sc), the
   mean over the top-5 by I_sc, and the lowest of all 25 decoys, plus the
   relaxed starting-model RMSD. Produces the data behind Figure 2 and Table S4.
2. `02_score_affinity_validation.ipynb` — Spearman correlation between Rosetta
   score metrics and experimental log(IC50)/log(Kd), pooled and per allele.
   Produces the data behind Figure 3 and Tables S3/S6.
3. `03_pymol_figures.ipynb` — PyMOL renders of representative decoy/crystal
   overlays (Figure 2A–C).
4. `04_figure1_panels.ipynb` — dataset composition panels (Figure 1).
5. `05_figure2_panels.ipynb` — assembles Figure 2 from the notebook-01 outputs.
6. `06_figure3_panels.ipynb` — assembles Figure 3 (I_sc) and Supplementary
   Figure S1 (reweighted_sc).
7. `07_supplemental_tables.ipynb` — builds Supplementary Tables S1–S6 and writes
   them to a single .xlsx.

## Score metrics

Three Rosetta score metrics are reported per decoy, each summarized as the best
(lowest) and mean across the 25-decoy ensemble:

- `I_sc` — interface score (peptide–MHC interaction energy). **Primary metric**;
  strongest association with experimental affinity.
- `reweighted_sc` — reweighted score used in prior applications of this pipeline.
- `total_score` — full-pose Rosetta energy (dominated by MHC internal energy).

All scores are in Rosetta Energy Units (REU); lower is more favorable.

## Helper scripts

- `parse_scorefiles.py` — parses per-pair `score.sc` files into a summary table
  of best/mean I_sc, reweighted_sc, total_score, and pep_sc.
- `rebuild_metadata.py` — merges the parsed score summaries into `metadata.csv`,
  adding the six score-summary columns.
- `attrition_counts.py` — reproduces the record-attrition funnel (Supplementary
  Table S2) from the raw IEDB download and the released metadata.
- `quick_metric_comparison.py` — pooled Spearman correlations for each score
  metric (Supplementary Table S6).
- `censoring_sensitivity.py`, `censoring_diagnostic.py` — robustness checks on
  the assay-detection-limit censoring rule.

## Paths

Data paths are hardcoded at the top of each notebook/script and reflect the
authors' cluster layout; edit them for your environment. Key locations:

- Released dataset: `.../IEDB_data_clean/huggingface/`
- Per-pair Rosetta outputs: `.../IEDB_data_clean/pdb/{allele}/{peptide}/`
  (25 decoys named `{peptide}_input_{0001..0025}.pdb` plus a `score.sc`
  scorefile; the `description` column matches each PDB filename stem)
- Template database: `.../MHC_database/templates/`

## Chain conventions

- **Modeled structures**: chain A = MHC binding cleft (α₁/α₂ domains,
  ~180 residues); chain B = peptide.
- **Experimental crystal structures**: chain IDs vary; the peptide chain is
  identified by sequence match (length 7–15 as fallback), and the MHC heavy
  chain as the remaining chain in the 170–290 residue range.

## Reproducing from the HuggingFace release

The per-decoy `pdb/` directory is not part of the public release. To reproduce
the structural analyses, extract per-decoy PDBs from the released silent files:

```bash
# from inside structures/{allele}/
extract_pdbs.linuxgccrelease -in:file:silent {peptide}.silent
```

## Dependencies

`pandas`, `numpy`, `scipy`, `matplotlib`, `biopython`, `tqdm`, `openpyxl`
(for the .xlsx supplementary tables). See `requirements.txt`. Rosetta is
required only for structure generation and `extract_pdbs`, not for the
analysis notebooks. `utils/` holds shared I/O, plotting, and structure helpers.
