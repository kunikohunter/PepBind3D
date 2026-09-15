# PepBind3D: curation and validation code

Code for the Scientific Data descriptor pairing curated IEDB peptide-HLA class I
binding affinities with Rosetta FlexPepDock structural ensembles: notebooks that
produce the manuscript figures and tables, the analysis scripts behind every
reported number.

Dataset: https://huggingface.co/datasets/kunikohunter/PepBind3D

**Current scope:** 95 HLA-A/B/C alleles, 112,378 peptide-allele pairs, 25 decoys
each (2,809,450 structures, ~124 GB), 118,751 measurement rows. Structures came
from two Rosetta builds, 3.14 for the initial batch, 2024.09+release.06b3cf8
for the alleles added later; `metadata.csv`'s `source_version` column says which.

## Structure generation: where to start reading

`IEDBTestPipeline_ACCRE.py` (~1,700 lines, 16 top-level functions, no classes) is
the production script. It is long, and most of it is not on the path that
produced this dataset, so start here rather than at the top of the file.

**Two independent jobs live in one file.** Curation and structure generation
share it but do not call each other:

| | entry point | what it does |
|---|---|---|
| curation | `get_peplist` → `clean_peplist` (+ `standardize_columns`, `rename_columns`, `flatten_columns`) | reads the IEDB bulk export, filters to quantitative IC50/KD on HLA-A/B/C, deduplicates per allele, writes the peptide lists |
| structure generation | `thread_all` → `thread_template` | builds one prepacked starting model per peptide, and optionally the SLURM + options files for the FlexPepDock production run |

**The generation path, in order** (all inside `thread_template`):

1. `HLA_db.MHCdatabase.get_peptide_template` picks the threading template:
   every same-length peptide from every allele of the same **gene**, ranked by
   BLOSUM62 similarity to the target. See the caveat below.
2. `SimpleThreadingMover` mounts the query peptide on the template peptide's
   backbone.
3. `DeleteRegionMover` trims the receptor to the α₁/α₂ cleft and removes
   template overhang.
4. `add_NCAA` substitutes non-canonical residues where a `.params` file is
   supplied.
5. the receptor is added and the complex relaxed, `FastRelax`, 5 repeats,
   ref2015.
6. FlexPepDock **prepack** produces the starting model that docking consumes.

Refinement itself is a separate SLURM array, not this script:
`-pep_refine -nstruct 25 -ex1 -ex2aro`.

**Reading the rest.** `postprocessing_affinity`, `build_scorefile` and
`stats_from_scorefile` are post-run bookkeeping over completed docking output.
`make_batch` / `create_batch_list` / `safe_thread_all` and the `--threads`,
`--slurm_setup`, `--batch_index` flags are the SLURM batch-array layer added for
this dataset; the generation logic underneath is unchanged from the published
version.

**Flags that change the output.** `--ignore_epitope_match` enables
self-template exclusion (`omit=["self"]`). It is **off** by default, which is how
the released structures were generated; `regeneration/` re-runs the
crystal-matched validation subset with it **on**. `--find_worst_template`
inverts the template ranking and exists for diagnostics only.

`HLA_db.py` holds the template database. `MHCdatabase` builds and queries it;
`residueSelect` and `chainSelect` are BioPython `PDBIO.Select` subclasses, so
their `accept_*` methods are called by BioPython rather than from this codebase.

## Provenance: the published prior version

The pipeline began as the code published with Bloodworth N, Chen W, Hunter K,
Patrick D, et al. *Posttranslationally modified self-peptides promote
hypertension in mouse models.* J Clin Invest. 2024;134(16):e174374.
doi:10.1172/JCI174374, https://github.com/meilerlab/discovery-self-peptides-hypertension
(`code/`).

This dataset used an adapted copy, included here because it, not the published
version, is what produced these structures. The adaptation adds SLURM
batch-array execution (`--batch_index`, `--slurm_setup`, `--threads`,
`create_batch_list`, `safe_thread_all`) and column-standardization helpers, and
points at a newer Rosetta tree. The generation logic, threading order,
template selection, docking protocol, is unchanged from the published version.

Structures came from two Rosetta builds: **3.14** for the initial batch of
alleles and **2024.09+release.06b3cf8** for those added subsequently. Both share
the FlexPepDock protocol and the ref2015 weights, and `metadata.csv`'s
`source_version` column records which produced each pair.

### Template selection, and the self-templating caveat

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

Note: notebooks 01 and 03 reconstruct the template choice with a helper that
filters to the same *allele*, which the real selector does not, their
`template_identity` values may not name the template actually used.

## Analysis scripts

Every number in the manuscript and response letter comes from one of these, and
each has a `--self-test` with an analytically known answer.

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

Defaults are the authors' layout, so the code reproduces the published analyses
unchanged when run in place. The **notebooks** still carry absolute paths in
their path constants and in their stored output: they are the record of how the
figures were produced rather than reusable tooling, so edit the constants at the
top of each if you re-run them.

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
