# Re-dock the 1,162 content-defective v1 pairs (ACCRE)

## Scope: v1 only, 1,162 pairs — definitive, not an estimate

Both halves of the release have now been screened with the same three-class
test (chains == {A,B} AND chain B sequence == the curated peptide):

| batch | screened | defective |
|---|---|---|
| v1 (Tungsten) | 49,268 | **1,162 (2.36%)** |
| v2 (ACCRE) | 63,110 | **0** |

v1 defects: **984 extra_chain**, **178 truncated_peptide**. Spread over 32 of
37 alleles, so not localised. 10-mers fail at 7.1% against 2.0% for 9-mers.

**799 of the 1,162 are invisible to any score-based check** — only 363 show the
all-zero interface-term signature. That is why this needed a content screen.

## Why v2 is clean, and what it implies for the re-dock

Almost certainly the selector patch: `HLA_db.py` was modified 11 May, between
the v1 and v2 production runs. Whatever the mechanism, the empirical fact is
what matters here — **the current pipeline produced 0 defects across 63,110
pairs**, so re-running these 1,162 through it should simply come out clean.

That is the prediction to check early rather than at the end: if staging starts
rejecting these pairs en masse, the current pipeline does *not* fix the
underlying problem and we need the root-cause fix below before burning the
hours. **Please run one allele first and confirm the staging rejection count is
~0 before submitting the full array.**

## Root cause (for your awareness, and in case the prediction fails)

`database.info` declares an epitope length, but the template *file* often has a
different number of resolved residues in that chain. Measured on 400 templates,
**63 (16%) disagree**. The selector filters candidates on the DECLARED length
(`Seq_Length == len(query)`), then threading maps the query onto the file's
actual residues:

- template file has MORE resolved residues than declared → surplus left behind
  as its own chain → `extra_chain` (984 cases)
- template file has FEWER → query truncated → `truncated_peptide` (178 cases)

Confirmed instance: **1JF1** is declared as the 10-mer `ELAGIGILTV` but has 11
resolved residues, and 1JF1 is exactly the template that produced the spurious
2-residue chain B on `A*02:01/EAAGIGILTV` in the 24-pair job.

The durable fix is to filter templates on the **resolved** residue count in the
template file rather than the declared `Seq_Length`. Only needed if the
prediction above fails.

## Settings — note these are the OPPOSITE of the 24-pair job

These are **release data**, so they must match the rest of the release:

- **NO `--ignore_epitope_match`.** Production default, self-matching templates
  allowed, exactly as the other 111,216 v1 pairs and all 63,110 v2 pairs were
  built. (The 24-pair validation job is the opposite case and needs the flag —
  do not let the two cross.)
- **Build 2024.09+release.06b3cf8**, same as v2 and the 201.
- Threading: `--IEDBquery skip --buildFasta <allele>_redock.fasta 0
  --setAllele "<allele>" --rosetta <ROSETTA_2024.09_MAIN> --threads <n>`
- Refinement: prepack, then `-pep_refine -nstruct 25 -ex1 -ex2aro`
- **Content gate at staging** (your `stage_redock.py`) — assert chains == {A,B}
  and chain B sequence == the pair's peptide. This is the thing that must not be
  skipped; it is the only reason the 201 and the 24 are clean by construction.

## Cost

1,162 pairs x 25 decoys = 29,050 decoys, 1 CPU per task, ~2.5 h walltime each.
~2,900 CPU-hours. At the 1,500-core concurrency available this is comfortably an
overnight run.

## Inputs

Generated on Tungsten, 18 KB total:

```
.../IEDB_data_clean/IEDB_validation/redock_1162/fastas/*_redock.fasta   # 32 files
.../IEDB_data_clean/IEDB_validation/redock_1162/allele_manifest_redock.csv
.../IEDB_data_clean/IEDB_validation/decoy_screen_v1/redock_1162_pairs.csv
```

`redock_1162_pairs.csv` carries `allele, allele_dir, peptide, verdict, detail`
so each pair's original defect is on record.

```bash
# from Tungsten (KH runs this; SSH from Tungsten needs interactive auth)
rsync -avP \
  $PEPBIND3D_DATA/IEDB_validation/redock_1162/ \
  $PEPBIND3D_DATA/IEDB_validation/decoy_screen_v1/redock_1162_pairs.csv \
  $CLUSTER_LOGIN:$PEPBIND3D_CLUSTER/main_project/data/redock_1162/
```

## Return and integration

Return to `.../IEDB_data_clean/IEDB_validation/redock_1162/pdb/{allele_dir}/{peptide}/`
with 25 decoys plus `score.sc` whose `description` matches the PDB stems.

On Tungsten these then go through `release/parse_scorefiles.py` and
`release/convert_to_silent.py` — the converter now has the same content gate, so
a defective pair cannot reach a silent file even if one slipped through
staging — and replace the defective silents in `release_v2_final/structures/`.
`analysis/screen_decoy_content.py` is re-run afterwards to confirm 0 defects.

## Effect on the manuscript

Small, and worth stating plainly: excluding all 1,162 moves score-affinity rho
by **+0.016 (IC50)** and **+0.002 (KD)**, and they are 0.98% of measurement
rows. Validation 1 is untouched, since it uses the regeneration re-dock rather
than released decoys. The 2 crystal-matched pairs among the 1,162 are
`A0201/CLGGLLTMV` and `A0201/TLACFVLAAV` — independently, the same 2 that failed
peptide correspondence when RMSD was computed on the released decoys.

Re-docking rather than dropping keeps the pair count at 112,378 (+183 orphans
= 112,561) rather than losing 1%.
