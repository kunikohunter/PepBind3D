# PepBind3D — curation and validation code

Code for the data descriptor *PepBind3D: Curated peptide-HLA class I binding
affinities with Rosetta structural ensembles*. This repository contains the
IEDB curation pipeline, the structure-generation scripts, and the notebooks
that produce the manuscript figures and supplementary tables.

The released dataset is on HuggingFace:
https://huggingface.co/datasets/kunikohunter/IEDB_pHLA_binding_data

## Curation and structure generation

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
