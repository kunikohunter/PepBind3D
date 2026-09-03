#!/usr/bin/env python3
"""
make_threading_inputs.py
=========================

Reconstruction of a lost pipeline stage: turns a peptide target list
(hlac_threading_targets.csv) into per-allele batch FASTA files plus a
threading_commands.sh command file, in the same format used by the
original (now-lost) generator that produced
/home/huntek1/main_project/scripts/slurm_batch_commands.sh.

The original generator script no longer exists. Every format decision
below was reverse-engineered from surviving artifacts it PRODUCED, or
from the consumer script that reads its output. Each decision cites the
specific surviving file/line it came from. Anything not directly provable
from a surviving artifact is called out explicitly in a "GUESS" comment.

-----------------------------------------------------------------------
FORMAT DECISIONS AND THEIR EVIDENCE
-----------------------------------------------------------------------

1. Command template
   Source: /home/huntek1/main_project/scripts/slurm_batch_commands.sh
            (2,009 lines, 37 alleles), e.g. line 1:
       python IEDBTestPipeline_ACCRE.py --IEDBquery skip --setAllele B*38:01 \
           --buildFasta ../data/IEDB_data/B*38:01/B*38:01_batch1.fasta 0
   Confirmed constant across the whole file (grep -c "HLA-" == 0; the
   trailing positional argument is always the literal "0" for every one
   of the 2,009 lines: awk '{print $NF}' | sort -u -> single value "0").
   The trailing "0" is args.buildFasta[1] in IEDBTestPipeline_ACCRE.py's
   fasta_rebuild() (lines 466-500): int(seq_toAdd)==0 short-circuits and
   just uses the fasta file as-is (no CSV re-merge). We reproduce it
   literally as "0" on every line, matching every observed line.

2. Batch numbering: unpadded, 1-indexed ("batch1", "batch2", ... "batchN")
   Source: slurm_batch_commands.sh itself -- grep -o 'batch[0-9]*\.fasta'
   shows values up to batch375 (A*02:01, 375 lines) with NO zero-padding
   anywhere in the file at any magnitude (batch1, batch41, batch375, etc).
   Also matches IEDBTestPipeline_ACCRE.py make_batch() (lines 580-629),
   which builds batch dirs/files with an unpadded f-string:
       batch_dir = os.path.join(output_dir, f'{fasta_base_name}_batch{b_dir_count}')
       tmp_fn = os.path.join(batch_dir, f'batch{b_dir_count}.fasta')
   NOTE (flag): a separate surviving artifact,
   .../full_peptide_docking/HLA-A0201/fasta_files/HLA-A0201_batch001.fasta,
   uses 3-digit zero-padding. This is inconsistent with both
   slurm_batch_commands.sh and make_batch()'s own f-string, and is judged
   to come from a different/later script variant. We follow the larger,
   internally-consistent, code-corroborated convention (unpadded) as the
   reconstruction target. This is the one place where surviving artifacts
   genuinely disagree -- flagged in the report, not silently resolved.

3. Batch size = 25
   Source: verified mathematically against slurm_batch_commands.sh line
   counts per allele vs. per-allele source FASTA record counts:
     - B*38:01: HLA-B_38_01_IEDB_data.fasta has 142 records (grep -c ">");
       ceil(142/25) = 6; slurm_batch_commands.sh has exactly 6 B*38:01
       lines (batch1..batch6).
     - B*39:01: HLA-B_39_01_IEDB_data.fasta has 1307 records;
       ceil(1307/25) = 53; slurm_batch_commands.sh has exactly 53
       B*39:01 lines.
   Also matches the code default: make_batch(fasta_fn, size=None) ->
       size = size or 25   # IEDBTestPipeline_ACCRE.py:583
   and create_batch_list() (line ~568) hard-codes make_batch(fasta_fn,
   size=25) as well. Default batch size below is therefore 25.

4. Allele string form: "<LETTER>*<NN>:<NN>" (e.g. "C*01:02"), IDENTICAL
   in --setAllele and inside the file path.
   Source: slurm_batch_commands.sh, e.g. line 100:
       --setAllele B*07:02 --buildFasta ../data/IEDB_data/B*07:02/B*07:02_batch41.fasta 0
   This exactly matches the "allele" column of the input CSV
   (hlac_threading_targets.csv, e.g. "C*01:02") with NO transformation
   needed -- the CSV already uses the same convention as --setAllele.
   (The separate HLA-<LETTER>_<NN>_<NN> directory-naming convention seen
   under /home/huntek1/Data/MHC_database/IEDB_data/ is a different,
   upstream naming scheme used only for the raw per-allele IEDB source
   dumps; it is not the convention used by slurm_batch_commands.sh or by
   our input CSV, so no mapping is needed here.)

5. FASTA header format for the batch files this stage produces:
   plain ">PEPTIDE" header, peptide sequence on the next line, NO
   "_mut_N" NCAA-variant records.
   Source: the per-allele source FASTA files that feed the ACCRE
   consumer (e.g. HLA-B_38_01_IEDB_data.fasta) use ">EpitopeID" as header
   with the raw peptide as the sequence line -- a bare identifier + bare
   sequence, no extra annotation. The two surviving downstream BATCH
   fasta files named in the task prompt (HLA-A0201_batch001.fasta,
   HLA-B2705/fasta_files/...) DO additionally show, after each canonical
   ">PEPTIDE" record, a ">PEPTIDE_mut_N" NCAA-substituted variant.

   CORRECTION (from independent review): an earlier version of this
   docstring wrongly attributed those "_mut_N" records to add_NCAA() in
   IEDBTestPipeline_ACCRE.py. That is false -- add_NCAA() rewrites an
   existing fasta in place and never emits a new "_mut_" header. The
   actual producer of "_mut_N" records, confirmed by reading it directly,
   is mutate_peptide()/edit_fasta() in
   /home/huntek1/main_project/scripts/thread_IEDB_peptides_cluster.py
   (lines 15-90). That script drives a DIFFERENT consumer script,
   IEDBTestPipeline.py (not IEDBTestPipeline_ACCRE.py), with extra
   --receptor/--params flags, for the full_peptide_docking project --
   i.e. exactly the "RELATED project" the task prompt told us to treat
   as "format reference, not gospel" for FASTA headers, not as the
   authority on the command template (that authority was explicitly
   assigned to slurm_batch_commands.sh / IEDBTestPipeline_ACCRE.py).

   DECISION (still a GUESS, flagged explicitly, and the single most
   consequential judgment call in this script): we emit ONLY the
   canonical ">PEPTIDE" record, no "_mut_N" variants, because our
   command template and consumer (IEDBTestPipeline_ACCRE.py via
   --buildFasta) is anchored to the ACCRE pipeline, which -- as far as
   the surviving code in that script shows -- has no mutate_peptide-style
   step of its own. hlac_threading_targets.csv has no Epitope_ID column
   either, so we could not reproduce the source fasta's own header
   convention exactly regardless.
   RISK: hlac_threading_targets.csv carries a "best_template_pdb" column,
   which mirrors thread_IEDB_peptides_cluster.py's "--receptor" argument
   much more closely than anything in the ACCRE flow. If the actual lost
   generator for THIS dataset was closer to thread_IEDB_peptides_cluster.py
   than to slurm_batch_commands.sh, then _mut_N variants SHOULD be
   included and the command template should carry --receptor/--params.
   We followed the task prompt's explicit authority ordering (A over B)
   rather than this column-name inference, but flag it here as a real,
   unresolved possibility that could make this script's output diverge
   from the true lost generator.

6. Deduplication: peptides are deduplicated per allele before batching.
   GUESS (flagged): the task instructions explicitly ask for this: no
   surviving artifact directly proves the original generator deduplicated,
   since the source IEDB fasta files are keyed by unique Epitope_ID and
   are not observed to contain literal duplicate peptide sequences within
   one allele. We deduplicate on (allele, peptide) as instructed, using
   first-seen order to keep batch numbering deterministic.

7. Command file path convention.
   slurm_batch_commands.sh uses paths relative to the scripts/ directory
   it was run from: "../data/IEDB_data/<ALLELE>/<ALLELE>_batchN.fasta",
   consistent with "python IEDBTestPipeline_ACCRE.py" also being resolved
   relative to that same scripts/ cwd -- i.e. every path on the original
   line resolves from one single, consistent working directory.
   Since this reconstruction writes to a brand-new location
   (/data/p_csb_meiler/huntek1/hlac_generator/output/), a relative path
   from "scripts/" no longer makes sense, and there is no single cwd from
   which BOTH a relative "IEDBTestPipeline_ACCRE.py" and a relative
   --buildFasta path under --outdir would resolve together. FIX (was a
   real bug in an earlier version of this script, caught by independent
   review: it emitted a path relative to --outdir while leaving the
   script invocation relative to scripts/, which cannot both resolve from
   one cwd): we use an ABSOLUTE path for --buildFasta instead. This
   preserves the command SHAPE exactly (flag order, allele string, "0"
   trailing arg) while being runnable from anywhere.

-----------------------------------------------------------------------
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict, defaultdict

DEFAULT_BATCH_SIZE = 25  # see decision 3 above


def load_targets(csv_path):
    """Read hlac_threading_targets.csv and group deduplicated peptides by allele.

    Columns: allele,peptide,best_template_pdb
    Returns: OrderedDict[allele] -> list of unique peptides, in first-seen order.
    """
    by_allele = OrderedDict()
    seen = defaultdict(set)
    total_rows = 0
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"allele", "peptide", "best_template_pdb"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Input CSV missing expected columns: {missing} "
                f"(found: {reader.fieldnames})"
            )
        for row in reader:
            total_rows += 1
            allele = row["allele"].strip()
            peptide = row["peptide"].strip()
            if not allele or not peptide:
                continue
            if allele not in by_allele:
                by_allele[allele] = []
            if peptide not in seen[allele]:
                seen[allele].add(peptide)
                by_allele[allele].append(peptide)
    return by_allele, total_rows


def batch_list(items, size):
    """Split items into consecutive chunks of `size`, 1-indexed batch numbers."""
    batches = []
    for i in range(0, len(items), size):
        batches.append(items[i : i + size])
    return batches


def write_fasta_batch(path, peptides):
    """Write a batch FASTA file: '>PEPTIDE' header, peptide as sequence line.
    See format decision 5 above."""
    with open(path, "w") as fh:
        for pep in peptides:
            fh.write(f">{pep}\n{pep}\n")


def build_command_line(allele, batch_num, fasta_rel_path):
    """Reproduce the exact command shape of slurm_batch_commands.sh.
    See format decisions 1 and 7 above."""
    return (
        f"python IEDBTestPipeline_ACCRE.py --IEDBquery skip "
        f"--setAllele {allele} --buildFasta {fasta_rel_path} 0"
    )


def run(csv_path, outdir, batch_size, dry_run):
    by_allele, total_rows = load_targets(csv_path)

    n_alleles = len(by_allele)
    n_unique_peptides = sum(len(v) for v in by_allele.values())
    n_batches_total = sum(
        len(batch_list(v, batch_size)) for v in by_allele.values()
    )

    print(f"Input CSV: {csv_path}")
    print(f"  data rows read:      {total_rows}")
    print(f"  unique alleles:      {n_alleles}")
    print(f"  unique peptides:     {n_unique_peptides} (deduplicated per allele)")
    print(f"  batch size:          {batch_size}")
    print(f"  total batch files:   {n_batches_total}")
    print(f"  total command lines: {n_batches_total}")

    if dry_run:
        print("\n--dry-run: no files written. Per-allele breakdown:")
        for allele, peptides in by_allele.items():
            nb = len(batch_list(peptides, batch_size))
            print(f"    {allele}: {len(peptides)} unique peptides -> {nb} batch(es)")
        return

    os.makedirs(outdir, exist_ok=True)
    iedb_data_dir = os.path.join(outdir, "IEDB_data")
    os.makedirs(iedb_data_dir, exist_ok=True)

    commands_path = os.path.join(outdir, "threading_commands.sh")
    written_fasta_files = 0
    written_peptides = 0

    with open(commands_path, "w") as cmd_fh:
        for allele, peptides in by_allele.items():
            allele_dir = os.path.join(iedb_data_dir, allele)
            os.makedirs(allele_dir, exist_ok=True)
            batches = batch_list(peptides, batch_size)
            for i, batch_peptides in enumerate(batches, start=1):
                fasta_name = f"{allele}_batch{i}.fasta"
                fasta_path = os.path.join(allele_dir, fasta_name)
                write_fasta_batch(fasta_path, batch_peptides)
                written_fasta_files += 1
                written_peptides += len(batch_peptides)

                # Absolute path: the emitted "python IEDBTestPipeline_ACCRE.py"
                # invocation is relative to wherever that script lives (per
                # slurm_batch_commands.sh, the scripts/ dir), while our
                # --buildFasta target lives under --outdir. There is no single
                # cwd from which both a relative script name and a relative
                # --buildFasta path resolve correctly, so we use an absolute
                # path here (an adaptation, not a claim about the original's
                # literal path text -- see format decision 7).
                fasta_abs_path = os.path.abspath(fasta_path)
                cmd_fh.write(
                    build_command_line(allele, i, fasta_abs_path) + "\n"
                )

    print(f"\nWrote {written_fasta_files} FASTA files "
          f"({written_peptides} total peptide records) under {iedb_data_dir}")
    print(f"Wrote {n_batches_total} command lines to {commands_path}")

    if written_peptides != n_unique_peptides:
        print(
            f"WARNING: written peptide count ({written_peptides}) != "
            f"unique peptide count ({n_unique_peptides})",
            file=sys.stderr,
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstructed generator: turns hlac_threading_targets.csv into "
            "per-allele batch FASTA files and a threading_commands.sh, "
            "matching the format of the (lost) generator behind "
            "slurm_batch_commands.sh. See module docstring for full "
            "format-derivation evidence and flagged assumptions."
        )
    )
    parser.add_argument(
        "--csv",
        default="/home/huntek1/main_project/scripts/hla_c/hlac_threading_targets.csv",
        help="Input target CSV (columns: allele,peptide,best_template_pdb).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory. Will contain IEDB_data/<ALLELE>/<ALLELE>_batchN.fasta "
             "and threading_commands.sh. Must be a NEW directory; nothing outside "
             "it is touched.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Peptides per batch FASTA file (default: {DEFAULT_BATCH_SIZE}, "
             f"reverse-engineered from make_batch()'s default and from "
             f"ceil(n_peptides/25) matching observed slurm_batch_commands.sh "
             f"batch counts per allele).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report peptide/batch/allele counts without writing any files.",
    )
    args = parser.parse_args()

    run(args.csv, args.outdir, args.batch_size, args.dry_run)


if __name__ == "__main__":
    main()
