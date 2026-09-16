---
license: cc-by-4.0
task_categories:
  - other
tags:
  - biology
  - protein-structure
  - MHC
  - HLA
  - peptide
  - binding-affinity
  - rosetta
  - immunology
pretty_name: 'PepBind3D: pHLA Binding Affinity with Rosetta FlexPepDock Structures'
size_categories:
  - 100K<n<1M
---

# PepBind3D

**Curated peptide-HLA class I binding affinities paired with Rosetta FlexPepDock
structural ensembles.** 112,561 peptide-allele pairs across 95 HLA-A, -B and -C
alleles, each with 25 docked decoy structures and per-decoy interface energies.

**Authors:** Kuniko Hunter, Rocco Moretti, Jens Meiler, David G. Harrison
Vanderbilt University Medical Center / Vanderbilt University
**License:** CC BY 4.0
**Metadata source:** [Immune Epitope Database (IEDB)](https://www.iedb.org)

---

## Dataset Summary

The interaction between peptide antigens and class I human leukocyte antigens
(HLA) is a central determinant of CD8+ T cell recognition, and a foundational
target for cancer immunotherapy, vaccine design and autoimmune disease research.
Computational prediction of peptide-HLA (pHLA) binding has progressed rapidly
with sequence-based machine learning, but structure-aware approaches remain
limited by the scarcity of paired experimental affinity data and
three-dimensional structural information.

This dataset pairs the two. Every peptide-allele pair carries an experimental
binding measurement curated from the IEDB and an ensemble of 25 Rosetta
FlexPepDock decoy structures with their energy terms, so a model can be trained
on measured affinity and modeled structure together. Each pair provides a
*distribution* over peptide conformations rather than a single pose, which a
co-folding model does not produce an equivalent of.

Pairs measured by both IC50 and Kd assays appear as separate rows, giving
118,985 measurement rows over 112,561 pairs.

Intended uses:

- Training and benchmarking machine learning models for pHLA binding affinity
  prediction, particularly models that consume explicit structural features
- Evaluating computational docking and scoring methods against experimental
  binding data
- Downstream structural analysis: anchor residue mapping, conformational
  diversity, interface comparison between binders and non-binders

---

## Quickstart

The structures total ~155 GB, so start with `metadata.csv` (40 MB) and pull only
the silent files you need.

```bash
pip install huggingface_hub pandas
```

**1. Get the metadata table.**

```python
import pandas as pd
from huggingface_hub import hf_hub_download

path = hf_hub_download("kunikohunter/PepBind3D", "metadata.csv",
                       repo_type="dataset")
df = pd.read_csv(path)
print(len(df), "measurements over", df[["allele", "peptide"]].drop_duplicates().shape[0], "pairs")
```

**2. Filter to usable measurements.** Roughly 44% of rows sit at or above an
assay detection ceiling and are not quantitative affinities; drop them before
fitting anything (see [Known Limitations](#known-limitations-and-caveats)).

```python
CEILINGS = {"IC50": ({20000, 50000, 70000}, 70000),
            "Kd":   ({5000, 10000, 20000},  20000)}

def quantitative(df, assay):
    exact, top = CEILINGS[assay]
    sub = df[(df.measurement_type == assay) & ~df.flagged]
    return sub[~(sub.measurement_value.isin(exact) | (sub.measurement_value >= top))]

ic50 = quantitative(df, "IC50")      # 18,113 rows
print(ic50[["allele", "peptide", "measurement_value", "I_sc_best"]].head())
```

**3. Download the structures for one pair.** `pdb_dir` gives the path for every
row.

```python
row = ic50.iloc[0]
silent = hf_hub_download("kunikohunter/PepBind3D", row.pdb_dir, repo_type="dataset")
```

To pull a whole allele instead, use `allow_patterns`:

```python
from huggingface_hub import snapshot_download
snapshot_download("kunikohunter/PepBind3D", repo_type="dataset",
                  allow_patterns="structures/A0201/*")
```

**4. Extract PDB coordinates** from a silent file with Rosetta:

```bash
extract_pdbs.linuxgccrelease -in:file:silent A0201/G/GILGFVFTL.silent
```

---

## Dataset Structure

```
PepBind3D/
├── README.md
├── metadata.csv                        # Master table - one row per IEDB measurement
└── structures/
    ├── A0101/
    │   ├── F/
    │   │   └── FHEFLSSKL.silent
    │   ├── G/
    │   │   └── GILGFVFTL.silent
    │   └── ...
    ├── A0201/
    └── ...
```

Each silent file holds the full 25-decoy ensemble for a single peptide-allele
pair, with Rosetta energy scores embedded alongside the coordinates. Allele
directories use filesystem-safe notation (`A0101/` for `HLA-A*01:01`); the
`allele` column in `metadata.csv` carries the standard notation (`A*01:01`).
Within each allele, files are grouped by the peptide's first residue, so
`A*02:01` / `GILGFVFTL` is at `structures/A0201/G/GILGFVFTL.silent`. The
`pdb_dir` column gives the full path for every row, so nothing needs to be
reconstructed by hand.

There is one silent file per unique peptide-allele pair (112,561 files); the
118,985 metadata rows exceed this because 6,423 pairs carry both an IC50 and a
Kd measurement and appear as two rows linked to the same silent file.

---

## Metadata Fields

| Column | Description |
|---|---|
| `allele_iedb` | Allele name in IEDB format (e.g. `HLA-A_01_01`) |
| `allele` | Standard allele notation (e.g. `A*01:01`) |
| `allele_compact` | Filesystem-safe allele form used in paths (e.g. `A0101`) |
| `peptide` | Peptide amino acid sequence |
| `peptide_length` | Length of peptide in residues |
| `measurement_type` | Assay type. Exact string values `IC50` or `Kd` (case-sensitive) |
| `measurement_value` | Quantitative binding measurement |
| `measurement_units` | Units of measurement (nM) |
| `assay_method` | Assay method as reported in IEDB |
| `assay_response` | Full assay response description from IEDB |
| `pubmed_id` | PubMed ID of the source publication |
| `parent_protein` | Source protein of the peptide |
| `protein_accession` | UniProt/GenBank accession of source protein |
| `source_organism` | Organism of origin |
| `source_version` | Curation batch the pair came from, `v1` or `v2` |
| `flagged` | Whether the entry was flagged during IEDB data cleaning |
| `self_templated` | Whether the pair's own crystal structure was available as a threading template. **Exclude these when measuring structure-prediction accuracy** |
| `has_structures` | Whether a structural ensemble is available for this pair |
| `num_pdbs` | Number of decoy structures in the silent file (25) |
| `I_sc_best` | Best (lowest) interface score across the 25-decoy ensemble (REU). **Recommended primary metric** - see below |
| `I_sc_mean` | Mean interface score across the ensemble (REU) |
| `reweighted_sc_best` / `reweighted_sc_mean` | Best and mean reweighted score (REU) |
| `total_score_best` / `total_score_mean` | Best and mean total_score (REU) |
| `pep_sc_best` / `pep_sc_mean` | Best and mean peptide score (REU) |
| `rosetta_best_score` / `rosetta_mean_score` | Aliases of `total_score_*`, retained for backward compatibility |
| `pdb_dir` | Relative path to the peptide silent file within the dataset |

**Which score to use.** Each score is summarized as the best (lowest) and mean
across the 25-decoy ensemble:

- `I_sc` (**interface score**) - the peptide-MHC interaction energy in
  isolation. This is the recommended primary metric; it showed the strongest
  association with experimental affinity.
- `reweighted_sc` (**reweighted score**) - upweights peptide-relevant energy
  terms but retains full-pose energy.
- `total_score` - the complete Rosetta pose energy. Because it is dominated by
  the internal energy of the MHC receptor, which varies little among peptides
  presented by the same allele, it is the least sensitive of the three for
  ranking binding.

All scores are in Rosetta Energy Units (REU); lower is more favorable. **REU is
not a binding free energy.** These values rank poses within a modeling
framework; they are not thermodynamic quantities and should not be read as
predicted affinities.

---

## Dataset Statistics

| Metric | Count |
|---|---|
| Peptide-allele pairs (silent files) | 112,561 |
| Total measurements (metadata rows) | 118,985 |
| Unique alleles | 95 (HLA-A, -B, -C) |
| Unique peptides | 25,622 |
| Peptide lengths | 7-15 residues |
| IC50 measurements | 21,234 |
| Kd measurements | 97,751 |
| Censored measurements | 52,462 (44%) |
| Flagged entries | 15 |
| Self-templated pairs | 370 (0.33%) |
| Total decoy structures | 2,814,025 |
| Approximate size | ~155 GB |

---

## IEDB Data Curation

Binding affinity data were retrieved from a local copy of the IEDB bulk download
(`mhc_ligand_full.csv`, accessed April 14, 2025) and processed with a custom
Python pipeline (`IEDBTestPipeline.py`).

### 1. Retrieval and Assay Filtering

The full IEDB MHC ligand table was read in chunks and filtered to entries
belonging to HLA-A, HLA-B or HLA-C alleles with quantitative binding
measurements. Rows were retained only if they reported one of two assay response
types:

- **Half maximal inhibitory concentration (IC50)**
- **Dissociation constant (Kd)**

Variant labels in the IEDB (`dissociation constant KD (~EC50)`,
`dissociation constant KD (~IC50)`) were normalized to a single canonical label
(`dissociation constant (KD)`) prior to filtering. Entries missing either a
quantitative measurement value or assay units were dropped.

### 2. Per-Allele Deduplication

Records were grouped by allele. Within each allele, duplicate entries for the
same epitope (identified by IEDB Epitope IRI) were resolved separately for IC50
and Kd measurements, using the following rules in order:

| Scenario | Action |
|---|---|
| Two entries, one lacks a PubMed ID | Retain the entry with a PubMed ID; drop the other |
| Two entries, both lack PubMed IDs, values differ by ≤10 nM | Retain one entry (first occurrence) |
| Two entries, both lack PubMed IDs, values differ by >10 nM | Drop both entries |
| Two entries, both have PubMed IDs, values differ by <10 nM | Retain one entry (first occurrence) |
| Two entries, both have PubMed IDs, values conflict (≥10 nM difference) | Flag both for manual review |
| More than two entries | Drop entries lacking PubMed IDs first; if >2 remain, retain the entry closest to the median measurement value |
| Any other ambiguous case | Flag for manual review |

Entries that could not be unambiguously resolved were written to a separate
`flagged_IEDB_data.csv` per allele and excluded from the cleaned dataset. These
are the `flagged = True` entries in `metadata.csv`.

### 3. Sequence-Level Filtering

Peptide sequences containing the `+` character (used by IEDB to denote
non-canonical or modified amino acids) were excluded from structure generation.
Retained peptides consist of the 20 standard amino acids, 7 to 15 residues long.

---

## Computational Methods

Structures were generated with **Rosetta FlexPepDock**, 25 decoys per
peptide-allele pair. Each peptide sequence was threaded onto a length-matched
template peptide from a local MHC template database (`SimpleThreadingMover`),
trimmed, given its receptor, relaxed (`FastRelax`, 5 repeats, ref2015),
prepacked, and refined by flexible peptide docking
(`-pep_refine -nstruct 25 -ex1 -ex2aro`). No score cutoff was applied: all 25
decoys are retained.

Where an allele had no experimental receptor structure, the α1/α2 domains were
modeled with AlphaFold2 from the IPD-IMGT/HLA protein alignment.

All structures were generated with Rosetta 2024.09+release.06b3cf8. The
`source_version` column records which curation batch a pair came from, v1 or v2,
not a difference in Rosetta build.

FlexPepDock refinement holds the MHC backbone fixed, so decoys within a pair
differ only in the peptide; MHC Cα coordinates are identical across an ensemble.

---

## Validation

Two analyses in the accompanying manuscript establish what the structures and
scores support.

- **Structural accuracy.** 76 pairs have a matching experimental crystal
  structure. Measured on ensembles re-docked with self-matching templates
  excluded, the best-scoring decoy reaches a median peptide-backbone RMSD of
  1.21 Å (IQR 0.93-1.68, 87% within 2 Å); the best decoy of each ensemble
  reaches 0.99 Å.
- **Score-affinity relationship.** Spearman ρ between best-decoy `I_sc` and log
  affinity is 0.31 (IC50, n = 18,113) and 0.18 (Kd, n = 48,395), censored values
  excluded. The relationship is real but modest: these scores are intended as
  input features for a downstream model, not as standalone affinity predictors.
- **Binder discrimination.** `I_sc` separates censored from quantitative
  measurements with AUROC 0.678 (IC50) and 0.639 (Kd), so the censored
  measurements carry information rather than being merely excluded.

---

## Known Limitations and Caveats

**Binding measurements**

- IC50 values pile up at 20,000, 50,000 and 70,000 nM, and Kd values at 5,000,
  10,000 and 20,000 nM, consistent with assay detection ceilings. Treat these as
  censored (≥ the reported value) rather than exact.
- **Values above the highest ceiling are censored too**, and testing only for the
  exact ceiling values keeps them: 1,366 IC50 rows sit above 70,000 nM and 8,651
  Kd rows above 20,000 nM, including a 1,000,000 nM placeholder used for
  peptides with no measurable binding. Censored rows are ~44% of the dataset
  overall, and half of the Kd rows.
- **Kd pools three different assays**, distinguishable through `assay_method`.
  The competitive radioligand subset is ~1% censored; the two fluorescence
  subsets are 71-74% censored and centered about one log unit stronger. Stratify
  on `assay_method` or model the censoring explicitly.
- IC50 and Kd are not directly comparable and have not been converted between
  each other. Both are retained, with `measurement_type` recording which.

**Structures**

- All structures are computationally generated, and the best-scoring decoy is
  not guaranteed to be the native-like conformation.
- RMSD values embedded in silent files are computed against the threading
  template, not against any experimental reference, and are not a measure of
  model accuracy.
- `self_templated` marks 370 pairs (0.33%) whose own crystal structure was in
  the threading template library, so their structures were built from real
  coordinates of that exact peptide. They are the most accurate structures here,
  and they must be excluded when measuring structure-prediction accuracy or the
  result is inflated.

**Coverage**

- Allele representation is uneven: `A*02:01` alone accounts for 9.2% of pairs,
  and seven alleles have fewer than five.
- HLA-C is included but sparsely sampled (1,861 pairs) and limited to 9- and
  10-mer peptides, because the template library holds only 9- and 10-mer HLA-C
  structures.
- HLA class II alleles are not included, and peptides carrying
  post-translational modifications or non-standard residues were excluded before
  structure generation.

---

## Citation

If you use this dataset, please cite:

> [Manuscript citation - to be added upon publication]

**Dataset DOI:** [https://doi.org/10.57967/hf/9669](https://doi.org/10.57967/hf/9669)

Experimental binding data are sourced from the
[Immune Epitope Database (IEDB)](https://www.iedb.org), also available under
CC BY 4.0. Please cite IEDB as well:

> Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB): 2018
> update. *Nucleic Acids Research*. 2019;47(D1):D339-D343.
> https://doi.org/10.1093/nar/gky1006

Template and reference structures come from the **RCSB PDB** (public domain,
CC0).

---

## Code

Curation, structure generation and every validation analysis:
https://github.com/kunikohunter/PepBind3D

---

## License

Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You
are free to share and adapt the material for any purpose, provided appropriate
credit is given.
