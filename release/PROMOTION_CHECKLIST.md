# Promotion checklist: staging to HuggingFace

**Decision point, not a script.** Promoting replaces a DOI'd public record, and
this step is reserved for KH. Nothing here runs automatically.

## What changes

| | live HF record | staging (`release_v2_final/`) |
|---|---|---|
| silents | 49,268 | 112,561 |
| measurement rows | 49,488 | 118,985 |
| alleles | 37 (A, B) | 95 (A, B, C) |
| unique peptides | 14,077 | 25,622 |
| metadata columns | 26 | 30 |
| content-defective pairs | 1,162 | 0 |
| last written | 21 July | current |

The live record also carries the Kd measurements the curation label typo
dropped, and 1,162 pairs whose structures do not contain the peptide they
claim. It is not a version anyone should be pointed at.

## Before promoting

- [ ] **Nothing cites the current DOI.** The manuscript is under review and the
      dataset has not been referenced elsewhere, so this is a version bump, not
      a correction. Confirm no preprint or talk has pointed at `10.57967/hf/9669`.
- [ ] **Re-run the audits** (all currently pass):
      - `python3 analysis/screen_decoy_content.py --out-dir <dir>`, expect 0 defects
      - referential integrity: every `pdb_dir` resolves, no orphan silents
      - `release/add_release_columns.py --self-test`
- [ ] **Confirm the file to ship is `release_v2_final/metadata.csv`**, 30
      columns, `assay_pdb_id` dropped, `pubmed_id` as integers.
      `metadata_curated.csv` is the 16-column curation input and must NOT ship.
- [ ] **Decide on `pep_sc`.** Currently shipped (4 score metrics, not 3). The
      dataset card documents it; make sure the manuscript's Data Records matches.

## Promoting

- [ ] Copy `release_v2_final/metadata.csv` into the HuggingFace working tree,
      and place structures at `structures/{allele}/{first residue}/{peptide}.silent`.
      **The two trees differ on purpose:** `release_v2_final/` is flat, because
      that is what the analysis scripts read; the staging tree is sharded,
      because the Hub caps a directory at 10,000 entries and the largest allele
      holds 10,376. `pdb_dir` in `metadata.csv` records the sharded path.
      **Do not `cp` over existing files**: the staging tree is
      hardlink-assembled and copying over a path writes through the shared
      inode into the source tree. Unlink first, or hardlink into a clean
      directory.
- [ ] Replace `README.md` with `release/DATASET_CARD.md`.
- [ ] Upload with the Hub API, not `git push`. Three limits apply and each one
      only appears after the previous is solved: 1000 API requests / 5 min
      (killed a `git push` at 30% after 3.5 h), 128 commits / hour (so commit
      per allele, not per 125 files), and the 10,000-entry directory cap.
      `hf upload-large-folder` shrinks its batch when rate-limited, which makes
      the commit ceiling worse rather than better.
- [ ] Confirm the new version resolves over HTTPS and downloads without
      authentication.

## After promoting

- [ ] **Record the new DOI.** It replaces `10.57967/hf/9669` in three places:
      Data Records, the data citation in the reference list, and the response
      file (`[CONFIRM]` markers).
- [ ] Re-run `revisions/apply_docx_edits.py`, counts and statistics refresh
      from the analysis outputs automatically; only the DOI is manual.
- [ ] **Mint the code DOI.** Zenodo → GitHub integration, enable
      `kunikohunter/PepBind3D`, cut a release, cite the *version* DOI in Code
      Availability. This makes the code citation permanent regardless of
      repository state, which is what Reviewer 1's "not accessible" comment was
      really about.
- [ ] Read the Methods line number off the line-numbered manuscript and fill the
      two `[LINE]` placeholders in E3.

## Not doing, and why

- **No correction notice.** Nothing was published against the old version.
- **No migration guide.** The old version should not be used; the new one is a
  superset with corrected structures.
- **Old version not deleted.** HuggingFace keeps history; leaving it costs
  nothing and preserving the record is better than hiding it.
