# PepBind3D: curation and validation code

Code for the Scientific Data descriptor pairing curated IEDB peptide-HLA class I
binding affinities with Rosetta FlexPepDock structural ensembles: notebooks that
produce the manuscript figures and tables, the analysis scripts behind every
reported number.

Dataset: https://huggingface.co/datasets/kunikohunter/PepBind3D

**Scope:** 95 HLA-A/B/C alleles, 112,561 peptide-allele pairs, 25 decoys each
(2,814,025 structures, ~124 GB), 118,985 measurement rows.

## Structure generation

`pipeline/IEDBTestPipeline_ACCRE.py` produced these structures. It is our
adapted copy of the pipeline published with Bloodworth N, Chen W, Hunter K,
Patrick D, et al. *Posttranslationally modified self-peptides promote
hypertension in mouse models.* J Clin Invest. 2024;134(16):e174374,
doi:10.1172/JCI174374. We added batch-array execution and helpers for the
current IEDB schema; threading order, template selection and docking protocol
are unchanged. It came from a general-purpose pipeline, so much of it is
unrelated to this dataset.

Two jobs share the file and never call each other:

| | entry point |
|---|---|
| curation | `get_peplist`, then `clean_peplist` |
| structure generation | `thread_all`, then `thread_template` |

Inside `thread_template`: pick a template (`HLA_db.get_peptide_template`),
mount the peptide on its backbone (`SimpleThreadingMover`), trim the receptor to
the α₁/α₂ cleft (`DeleteRegionMover`), substitute non-canonical residues
(`add_NCAA`), add the receptor and relax (`FastRelax`, 5 repeats, ref2015), then
prepack. Refinement is a separate SLURM array: `-pep_refine -nstruct 25 -ex1
-ex2aro`.

`--ignore_epitope_match` excludes a template whose peptide matches the target.
It is off by default, which is how the released structures were built;
`regeneration/` re-runs the crystal-matched subset with it on.

`HLA_db.py` builds and queries the template database.

## Template selection, and the self-templating caveat

`HLA_db.MHCdatabase.get_peptide_template` pools every same-length peptide from
every allele of the **same gene** (not the same allele) and ranks them by
BLOSUM62 similarity to the target, highest first. Sequence identity and
resolution play no part in the ranking.

Because an *identical* peptide scores highest, a pair whose own crystal is in
the template database gets threaded onto itself. Self-exclusion exists
(`omit=["self"]`, matched on identical peptide sequence) but runs only under
`--ignore_epitope_match`, which the production runs did not pass. This affects
~0.1% of the release, but that 0.1% is exactly the set with crystal structures
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

## Analysis scripts

Every number in the manuscript and response letter comes from one of these.
Each takes a `--self-test` flag.

| | |
|---|---|
| `crystal_match.py` | matches release pairs to crystal structures in the template DB (76 pairs on the merged release; reproduces 52 on the earlier batch) |
| `crystal_rmsd.py` | RMSD for those pairs from the release silents, **upper bounds only**, see the self-templating caveat |
| `screen_decoy_content.py` | audits a decoy tree for content defects (spurious extra chain, truncated peptide) that no score-based check can see |
| `figure2g_template_identity.py` | rebuilds Figure 2G from the templates actually recorded in the threading logs |
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
| `release/add_release_columns.py` | adds `flagged`, `self_templated`, `has_structures`, `num_pdbs`, `pdb_dir` and the score columns to the merged metadata |
| `release/decoy_content.py` | the decoy content test, shared by the converter's gate and the audit script so the two cannot drift |
| `release/rebuild_metadata.py` | merges score summaries into `metadata.csv` |

## Score metrics

Per decoy, summarized as best (lowest) and mean over the 25-decoy ensemble, in
Rosetta Energy Units; lower is more favorable.

- `I_sc`, interface score. **Primary metric**, strongest association with affinity.
- `reweighted_sc`, score used in prior applications of this pipeline.
- `total_score`, full-pose energy, dominated by MHC internal energy.
- `pep_sc`, peptide-only score.

## Conventions

- **Alleles:** `A*02:01` in prose and metadata; `A0201` in filesystem paths.
- **RMSD:** peptide backbone (N, Cα, C, O) after superposition on the first 180
  MHC Cα atoms (α₁/α₂ domains). Always state this when reporting an RMSD.
- **Chains:** modeled structures are chain A = MHC cleft, chain B = peptide.
  Crystal chain IDs vary, the peptide chain is found by sequence match, the
  heavy chain as the remaining chain of 170–290 residues.
- **Affinity correlations:** Spearman on log₁₀ values with censored measurements
  excluded (IC50 20,000/50,000/70,000 nM; KD 5,000/10,000/20,000 nM).
- **KD measurements pool three IEDB assay labels** corresponding to different
  assays, distinguishable via `assay_method`. The competitive radioligand subset
  is ~1% censored; the two fluorescence subsets are 71–74% censored and centred
  about one log unit stronger. Stratify or model the censoring.

## Paths and dependencies

All filesystem roots live in `paths.py` and are overridable by environment
variable, so running this elsewhere means setting two variables rather than
editing every script:

```bash
export PEPBIND3D_DATA=/your/path/to/IEDB_data_clean    # measurements, structures, analysis outputs
export PEPBIND3D_MHC_DB=/your/path/to/MHC_database     # templates, database.info, raw IEDB export
export PEPBIND3D_ROSETTA=/your/rosetta/main            # only for structure generation
python3 paths.py                                       # prints the resolved roots and whether they exist
```

The defaults point at our own layout, so running the code in place reproduces
the published analyses unchanged. The notebooks set their paths in the first
cell; edit those if you re-run them.

Per-pair Rosetta outputs live at `$PEPBIND3D_DATA/pdb/{allele}/{peptide}/`
(`{peptide}_input_{0001..0025}.pdb` plus `score.sc`, whose `description` column
matches each PDB stem).

The per-decoy `pdb/` tree is not part of the public release. To reproduce the
structural analyses from it, extract from the released silent files:

```bash
# from inside structures/{allele}/
extract_pdbs.linuxgccrelease -in:file:silent {peptide}.silent
```

`pandas`, `numpy`, `scipy`, `matplotlib`, `biopython`, `tqdm`, `openpyxl`; see
`requirements.txt`. Rosetta is needed only for structure generation and
`extract_pdbs`, PyRosetta only for silent-file conversion, not for the analysis
notebooks. `utils/` holds shared I/O, plotting and structure helpers.
