# PepBind3D: curation and validation code

Code for the Scientific Data descriptor pairing curated IEDB peptide-HLA class I
binding affinities with Rosetta FlexPepDock structural ensembles: notebooks that
produce the manuscript figures and tables, and the analysis scripts behind
reported numbers.

Dataset: https://huggingface.co/datasets/kunikohunter/PepBind3D

**Scope:** 95 HLA-A/B/C alleles, 112,561 peptide-allele pairs, 25 decoys each
(2,814,025 structures, ~155 GB), 118,985 measurement rows.

## Structure generation

`pipeline/IEDBTestPipeline_ACCRE.py` produced these structures. It is an
adapted copy of the pipeline published with Bloodworth N, Chen W, Hunter K,
Patrick D, Palubinsky A, Phillips E, Roeth D, Kalkum M, Mallal S, Davies S,
Ao M, Moretti R, Meiler J, Harrison DG. *Posttranslationally modified
self-peptides promote hypertension in mouse models.* J Clin Invest.
2024;134(16):e174374, doi:10.1172/JCI174374. We added batch-array execution and helpers for the
current IEDB schema. Threading order, template selection and the docking
protocol are unchanged.

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

**Notebook 03 needs its own environment.** It is the only one that imports
`py3Dmol` and shells out to `pymol`, so run it in a PyMOL environment and the
rest in the analysis environment. The two are not interchangeable in either
direction: a PyMOL environment typically lacks `seaborn`, `openpyxl` and
`adjustText`, which notebooks 02, 06 and 07 need.

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
| `ensemble_diversity.py` | decoy-to-decoy vs decoy-to-crystal spread for the validation pairs |
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

## Score metrics

Per decoy, summarized as best (lowest) and mean over the 25-decoy ensemble, in
Rosetta Energy Units.

- `I_sc`, interface score. **Primary metric**, strongest association with affinity.
- `reweighted_sc`, score used in prior applications of this pipeline.
- `total_score`, full-pose energy, dominated by MHC internal energy.
- `pep_sc`, peptide-only score.

## Conventions

- **Alleles:** `A*02:01` in prose and metadata; `A0201` in filesystem paths.
- **RMSD:** peptide backbone (N, Cα, C, O) after superposition on the MHC Cα
  atoms, paired by global sequence alignment over the whole shared MHC
  region. The modeled receptor is trimmed to the α₁/α₂ cleft, so this is
  ~180 residues in practice, but the code does not impose a fixed window.
  Always state this when reporting an RMSD.
- **Chains:** modeled structures are chain A = MHC cleft, chain B = peptide.
  Crystal chain IDs vary, the peptide chain is found by sequence match, the
  heavy chain as the remaining chain of 170–290 residues.
- **Affinity correlations:** Spearman on log₁₀ values with censored measurements
  excluded (IC50 20,000/50,000/70,000 nM; KD 5,000/10,000/20,000 nM).
- **KD measurements pool three IEDB assay labels** corresponding to different
  assays, distinguishable via `assay_method`. They differ in location and in how
  much of each is censored, so stratify on `assay_method` rather than pooling
  blind. `analysis/kd_label_pooling.py` has the per-label numbers.

## Paths and dependencies

All filesystem roots live in `paths.py`. Set them for your machine:

```bash
export PEPBIND3D_DATA=/your/path/to/pepbind3d_data     # measurements, structures, analysis outputs
export PEPBIND3D_MHC_DB=/your/path/to/mhc_database     # templates, database.info, raw IEDB export
export PEPBIND3D_ROSETTA=/your/rosetta/main            # only for structure generation
python3 paths.py                                       # prints the resolved roots and whether they exist
```

The notebooks import the same roots in their first cell, so setting the
variables covers them too.

Notebook outputs are stored inside the `.ipynb`, so a cell that prints a path
writes it into the file and the next commit publishes it. Install the hook once
and that stops being something to remember:

```bash
ln -sf ../../tools/pre-commit .git/hooks/pre-commit
```

It rewrites paths in outputs to `<DATA_ROOT>`-style placeholders and re-stages
the notebook. A path hardcoded in a *code* cell fails the commit instead of
being rewritten, because the fix there is to import the root from `paths.py`.

Per-pair Rosetta outputs live at `$PEPBIND3D_DATA/pdb/{allele}/{peptide}/`
(`{peptide}_input_{0001..0025}.pdb` plus `score.sc`, whose `description` column
matches each PDB stem).

The per-decoy `pdb/` tree is not part of the public release. To reproduce the
structural analyses from it, extract from the released silent files:

```bash
# released silents are at structures/{allele}/{first residue}/{peptide}.silent
extract_pdbs.linuxgccrelease -in:file:silent A0201/G/GILGFVFTL.silent
```

See `requirements.txt`. Rosetta is needed only for structure generation and
`extract_pdbs`, PyRosetta only for silent-file conversion, neither for the
analysis notebooks. `utils/` holds shared I/O, plotting and structure helpers.
