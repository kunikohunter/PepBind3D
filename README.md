# PepBind3D — curation and validation code

Code for the Scientific Data descriptor pairing curated IEDB peptide-HLA class I
binding affinities with Rosetta FlexPepDock structural ensembles: notebooks that
produce the manuscript figures and tables, the analysis scripts behind every
reported number.

Dataset: https://huggingface.co/datasets/kunikohunter/PepBind3D

**Current scope:** 95 HLA-A/B/C alleles, 112,378 peptide-allele pairs, 25 decoys
each (2,809,450 structures, ~124 GB), 118,751 measurement rows. Structures came
from two Rosetta builds — 3.14 for the initial batch, 2024.09+release.06b3cf8
for the alleles added later; `metadata.csv`'s `source_version` column says which.

## Structure generation is not in this repository

It was published with Bloodworth N, Chen W, Hunter K, Patrick D, et al.
*Posttranslationally modified self-peptides promote hypertension in mouse
models.* J Clin Invest. 2024;134(16):e174374. doi:10.1172/JCI174374 — code at
https://github.com/meilerlab/discovery-self-peptides-hypertension (`code/`).

This dataset used an adapted copy. The adaptation adds SLURM batch-array
execution (`--batch_index`, `--slurm_setup`, `--threads`, `create_batch_list`,
`safe_thread_all`) and column-standardization helpers, and points at a newer
Rosetta tree. The generation logic — threading order, template selection,
docking protocol — is unchanged.

- `IEDBTestPipeline.py` / `IEDBTestPipeline_ACCRE.py` — curation and threaded
  starting models. Threading is SimpleThreadingMover → NCAA substitution → trim
  → add receptor → FastRelax (5 repeats, ref2015) → FlexPepDock prepack.
  Refinement is `-pep_refine -nstruct 25 -ex1 -ex2aro`.
- `HLA_db.py` — builds and queries the local MHC template database.

### Template selection, and the self-templating caveat

`HLA_db.MHCdatabase.get_peptide_template` pools every same-length peptide from
every allele of the **same gene** (not the same allele) and ranks them by
BLOSUM62 similarity to the target, highest first. Sequence identity and
resolution play no part in the ranking.

Because an *identical* peptide scores highest, a pair whose own crystal is in
the template database gets threaded onto itself. Self-exclusion exists
(`omit=["self"]`, matched on identical peptide sequence) but runs only under
`--ignore_epitope_match`, which the production runs did not pass. This affects
~0.1% of the release — but that 0.1% is exactly the set with crystal structures
to validate against, so:

- `regeneration/` re-runs the crystal-matched validation pairs with
  `--ignore_epitope_match`. **Validation 1 numbers must come from there, not
  from the release.** Measuring the released ensembles instead is optimistic by
  ~0.24 Å.

## Notebooks

Run in order; each reads `metadata.csv` and, where noted, per-pair score files.

| | |
|---|---|
| `01_structural_validation.ipynb` | peptide-backbone RMSD of decoys against matched crystals (Figure 2, Table S4) |
| `02_score_affinity_validation.ipynb` | Spearman correlation of score metrics against log affinity, pooled and per allele (Figure 3, Tables S3/S6) |
| `03_pymol_figures.ipynb` | decoy/crystal overlay renders (Figure 2A–C) |
| `04_figure1_panels.ipynb` | dataset composition (Figure 1) |
| `05_figure2_panels.ipynb`, `06_figure3_panels.ipynb` | figure assembly |
| `07_supplemental_tables.ipynb` | Supplementary Tables S1–S6 to one .xlsx |

Note: notebooks 01 and 03 reconstruct the template choice with a helper that
filters to the same *allele*, which the real selector does not — their
`template_identity` values may not name the template actually used.

## Analysis scripts

Every number in the manuscript and response letter comes from one of these, and
each has a `--self-test` with an analytically known answer.

| | |
|---|---|
| `crystal_match.py` | matches release pairs to crystal structures in the template DB (76 pairs on the merged release; reproduces 52 on the earlier batch) |
| `crystal_rmsd.py` | RMSD for those pairs from the release silents — **upper bounds only**, see the self-templating caveat |
| `ensemble_diversity.py` | decoy-to-decoy vs decoy-to-crystal spread, for the validation pairs |
| `recompute_affinity_112k.py` | censoring AUROC, pooled and per-allele Spearman, composition table |
| `kd_label_pooling.py` | recovers the original IEDB assay-response label for every KD row and tests whether the three are poolable |
| `affinity_baseline.py` | minimal sequence/structure affinity predictors |
| `attrition_counts.py` | record-attrition funnel (Table S2) |
| `censoring_sensitivity.py`, `censoring_diagnostic.py` | robustness of the censoring rule |
| `figure1_composition.py`, `figure2_diversity_panel.py` | figures regenerated on the merged release |

## Release helpers

| | |
|---|---|
| `release/parse_scorefiles.py` | per-pair `score.sc` → summary of best/mean `I_sc`, `reweighted_sc`, `total_score`, `pep_sc` |
| `release/convert_to_silent.py` | per-pair PDBs → one silent file, scores attached; hard-fails if the decoy count disagrees |
| `release/add_release_columns.py` | adds `flagged`, `has_structures`, `num_pdbs`, `pdb_dir` and the score columns to the merged metadata |
| `release/rebuild_metadata.py` | merges score summaries into `metadata.csv` |

## Score metrics

Per decoy, summarized as best (lowest) and mean over the 25-decoy ensemble, in
Rosetta Energy Units; lower is more favorable.

- `I_sc` — interface score. **Primary metric**, strongest association with affinity.
- `reweighted_sc` — score used in prior applications of this pipeline.
- `total_score` — full-pose energy, dominated by MHC internal energy.
- `pep_sc` — peptide-only score.

## Conventions

- **Alleles:** `A*02:01` in prose and metadata; `A0201` in filesystem paths.
- **RMSD:** peptide backbone (N, Cα, C, O) after superposition on the first 180
  MHC Cα atoms (α₁/α₂ domains). Always state this when reporting an RMSD.
- **Chains:** modeled structures are chain A = MHC cleft, chain B = peptide.
  Crystal chain IDs vary — the peptide chain is found by sequence match, the
  heavy chain as the remaining chain of 170–290 residues.
- **Affinity correlations:** Spearman on log₁₀ values with censored measurements
  excluded (IC50 20,000/50,000/70,000 nM; KD 5,000/10,000/20,000 nM).
- **KD measurements pool three IEDB assay labels** corresponding to different
  assays, distinguishable via `assay_method`. The competitive radioligand subset
  is ~1% censored; the two fluorescence subsets are 71–74% censored and centred
  about one log unit stronger. Stratify or model the censoring.

## Paths and dependencies

Paths are hardcoded at the top of each notebook/script for the authors' cluster
layout; edit for your environment. Per-pair Rosetta outputs live at
`.../IEDB_data_clean/pdb/{allele}/{peptide}/` (`{peptide}_input_{0001..0025}.pdb`
plus `score.sc`, whose `description` column matches each PDB stem); the template
database at `.../MHC_database/`.

The per-decoy `pdb/` tree is not part of the public release. To reproduce the
structural analyses from it, extract from the released silent files:

```bash
# from inside structures/{allele}/
extract_pdbs.linuxgccrelease -in:file:silent {peptide}.silent
```

`pandas`, `numpy`, `scipy`, `matplotlib`, `biopython`, `tqdm`, `openpyxl`; see
`requirements.txt`. Rosetta is needed only for structure generation and
`extract_pdbs`, PyRosetta only for silent-file conversion — not for the analysis
notebooks. `utils/` holds shared I/O, plotting and structure helpers.
