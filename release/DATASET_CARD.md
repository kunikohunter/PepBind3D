---
license: cc-by-4.0
language:
  - en
tags:
  - biology
  - immunology
  - protein-structure
  - peptide-mhc
  - binding-affinity
  - rosetta
pretty_name: Curated peptide-HLA class I binding affinities with Rosetta FlexPepDock structural ensembles
size_categories:
  - 100K<n<1M
---

# Curated peptide-HLA class I binding affinities with Rosetta FlexPepDock structural ensembles

Experimental peptide-HLA class I binding measurements curated from the Immune
Epitope Database, each paired with a 25-member Rosetta FlexPepDock structural
ensemble and per-decoy interface energy terms.

| | |
|---|---|
| peptide-allele pairs | 112,561 |
| measurements | 118,985 |
| HLA alleles | 95 (HLA-A, -B, -C) |
| unique peptides | 25,622 |
| decoy structures | 2,814,025 |
| size | ~155 GB |
| license | CC BY 4.0 |

## What this is for

Training and benchmarking structure-aware peptide-HLA binding predictors. Each
pair provides a *distribution* over peptide conformations rather than a single
pose, with physically-motivated interface energies attached to every decoy:
features a co-folding model does not produce an equivalent of.

## Layout

```
metadata.csv                                     one row per measurement
structures/{allele}/{P}/{peptide}.silent         25 decoys, Rosetta binary silent
```

`{P}` is the peptide's first residue, so `A*02:01` / `GILGFVFTL` is at
`structures/A0201/G/GILGFVFTL.silent`. The extra level exists because the Hub
caps a directory at 10,000 entries and the largest allele holds more than that;
sharding every allele by the same rule keeps the path predictable rather than
making one allele an exception.

Alleles use the filesystem-safe form in paths (`A0201`) and the standard form in
metadata (`A*02:01`), so the path for any row is:

```python
f"structures/{row.allele_compact}/{row.peptide[0]}/{row.peptide}.silent"
```

Extract individual PDBs with Rosetta's `extract_pdbs`:

```bash
extract_pdbs.linuxgccrelease -in:file:silent A0201/G/GILGFVFTL.silent
```

## metadata.csv columns

**Identity** `allele_iedb`, `allele`, `allele_compact`, `peptide`,
`peptide_length`

**Measurement** `measurement_type` (IC50 or Kd), `measurement_value`,
`measurement_units` (nM), `assay_method`, `assay_response`, `pubmed_id`

**Provenance** `parent_protein`, `protein_accession`, `source_organism`,
`source_version`, `flagged`, `self_templated`

**Structures** `has_structures`, `num_pdbs`, `pdb_dir`

**Scores**, best (lowest) and mean over the 25 decoys, in Rosetta Energy Units:
`I_sc_*` (interface), `reweighted_sc_*`, `total_score_*`, `pep_sc_*`, plus
`rosetta_best_score` / `rosetta_mean_score` as legacy aliases of `total_score_*`.

## Using it correctly

**REU is not a binding free energy.** The scores rank poses within a modeling
framework; they are not thermodynamic quantities and should not be read as
predicted affinities.

**Censored measurements.** IC50 at 20,000 / 50,000 / 70,000 nM and Kd at
5,000 / 10,000 / 20,000 nM are assay detection ceilings, not measurements.
Treat them as "≥ this value" or exclude them. They are roughly 40% of rows.

Exclude **anything at or above the highest ceiling** as well, not only exact
matches to the three values: 70,000 nM for IC50, 20,000 nM for Kd. A further
8,651 Kd rows (8.9%) and 1,366 IC50 rows sit above their top ceiling at values
like 77,900 nM, and testing only for the three listed numbers keeps them.

**Kd pools three different assays.** The three IEDB dissociation-constant
labels correspond to different experimental readouts and are distinguishable
through `assay_method`. The competitive radioligand subset is ~1% censored; the
two fluorescence subsets are 71–74% censored and centered about one log unit
stronger. Stratify on `assay_method` or model the censoring explicitly.

**`self_templated` marks 370 pairs (0.33%) whose own crystal structure was in
the threading template library**, so their structures were built from real
coordinates of that exact peptide. They are the most accurate structures here,
and they must be excluded when measuring structure-prediction accuracy, or the
result is inflated.

**Allele coverage is uneven.** HLA-A\*02:01 alone is 9.2% of pairs; seven
alleles have fewer than five. HLA-C is limited to 9- and 10-mer peptides,
because the template library holds only 9- and 10-mer HLA-C structures.

**The receptor is rigid.** FlexPepDock refinement holds the MHC backbone fixed,
so decoys within a pair differ only in the peptide; MHC Cα coordinates are
identical across an ensemble.

## How it was made

Measurements come from the IEDB MHC ligand bulk export (`mhc_ligand_full.csv`,
downloaded 14 April 2025), filtered to quantitative IC50/Kd on HLA-A, -B and -C
and deduplicated per allele. Structures were built by threading each peptide
onto a length-matched template peptide from a local MHC template database
(`SimpleThreadingMover`), trimming, adding the receptor, relaxing (`FastRelax`,
5 repeats, ref2015), prepacking, and refining with FlexPepDock
(`-pep_refine -nstruct 25 -ex1 -ex2aro`). No score cutoff was applied: all 25
decoys are retained.

All structures were generated with Rosetta 2024.09+release.06b3cf8. The
`source_version` column records which curation batch each pair came from, v1 or
v2, not a difference in Rosetta build.

## Validation

- **Structural accuracy.** 76 pairs have a matching experimental crystal
  structure. Measured on ensembles re-docked with self-matching templates
  excluded, the best-scoring decoy reaches a median peptide-backbone RMSD of
  1.21 Å (IQR 0.93–1.68, 87% within 2 Å); the best decoy of each ensemble
  reaches 0.99 Å.
- **Score-affinity relationship.** Spearman ρ between best-decoy `I_sc` and
  log affinity is 0.315 (IC50, n = 18,113) and 0.179 (Kd, n = 48,395),
  censored values excluded.
- **Binder discrimination.** `I_sc` separates censored from quantitative
  measurements with AUROC 0.678 (IC50) and 0.639 (Kd).

## Sources and license

Binding measurements: **IEDB**, https://www.iedb.org, CC BY 4.0
(https://www.iedb.org/citation_v3.php). Every measurement retains its PubMed
identifier where IEDB provides one, so values remain attributed to the
publishing authors.

Template and reference structures: **RCSB PDB**, public domain (CC0).

This dataset is released under **CC BY 4.0**. Please cite both this dataset and
the IEDB.

## Code

Curation, structure generation and every validation analysis:
https://github.com/kunikohunter/PepBind3D
