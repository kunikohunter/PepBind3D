#!/usr/bin/env python3
"""
run_chai1_msa.py

MSA-ON wrapper around Chai-1 (chai_lab==0.6.1) for the pMHC structure
prediction benchmark. This is a SECOND, SEPARATE Chai-1 wrapper alongside
(not replacing) benchmark/chai1/run_chai1_leakage_safe.py.

WHY MSAs ARE ON BUT TEMPLATES STAY OFF
----------------------------------------
MSAs are sequence-derived: an MSA hit tells the model "sequences like this
exist," which is co-evolutionary/statistical signal, not 3-D coordinates of
any specific structure. It cannot hand the model an already-solved answer
structure. Templates are the opposite: a template hit hands the model an
actual 3-D structure (coordinates) to copy from, which -- for a benchmark
whose holdout is "structures released after 2024-01-01" -- is a direct
leak vector if that template happens to be (or closely match) the answer
structure itself, pulled live from present-day RCSB.

This is why this benchmark turns MSAs ON everywhere (its primary condition)
while keeping templates OFF/date-restricted everywhere. This wrapper
enables `--use-msa-server` (remote MMSeqs2/ColabFold MSA generation) but
explicitly HARD-DISABLES `--use-templates-server` and never sets
`--template-hits-path`, mirroring the leakage-safe wrapper's discipline for
the one channel that actually matters for leakage.

Verified in chai_lab 0.6.1 (`chai-lab fold --help`):
  --use-msa-server / --no-use-msa-server        [default: no-use-msa-server]
  --msa-server-url TEXT   [default: https://api.colabfold.com]
  --use-templates-server / --no-use-templates-server
                                                 [default: no-use-templates-server]
  --template-hits-path PATH

USAGE
-----
    python run_chai1_msa.py \\
        --peptide-seq SIINFEKL \\
        --mhc-seq SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRT \\
        --output-dir /data/p_csb_meiler/huntek1/benchmark/wrappers/test_chai1_msa/output

Or call `run_chai1_msa(peptide_seq, mhc_seq, output_dir)` directly from
Python.

REACHING 25 SAMPLES
--------------------
chai-lab's own multi-sample mechanism is two nested knobs:
  --num-trunk-samples INTEGER   [default: 1]  -- independent trunk-module
                                                  samples (distinct
                                                  structure "hypotheses")
  --num-diffn-samples INTEGER   [default: 5]  -- diffusion samples drawn
                                                  PER trunk sample
Total structures emitted = num_trunk_samples * num_diffn_samples. With the
CLI defaults (trunk=1, diffn=5) this is exactly the 5-per-run observed in
earlier testing. To reach 25 per target, this wrapper sets BOTH to 5
(5 trunk samples x 5 diffusion samples/trunk = 25 total structures), rather
than only raising --num-diffn-samples, so the 25 samples also cover 5
independent trunk hypotheses instead of 5 diffusion draws off a single
trunk embedding -- closer in spirit to AF2's 5-models x 5-seeds and
Boltz's 25 independent diffusion samples.

NOTE: Chai-1 inference must be run on a GPU compute node (via sbatch/salloc),
never on the login node. This wrapper only builds the input FASTA and shells
out to `chai-lab fold`; it does not submit any Slurm job itself.

Requires the chai-lab conda env to be active, e.g.:
    source /data/p_csb_meiler/apps/miniforge3/bin/activate chai-lab
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

CHAI_LAB_BIN = "chai-lab"  # resolved via $PATH once the chai-lab env is active

# Same FASTA header convention as run_chai1_leakage_safe.py.
MHC_CHAIN_NAME = "mhc"
PEPTIDE_CHAIN_NAME = "peptide"


def build_fasta(peptide_seq: str, mhc_seq: str, output_dir: Path) -> Path:
    """Write the two-chain pMHC FASTA input Chai-1 expects."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "input.fasta"
    with open(fasta_path, "w") as fh:
        fh.write(f">protein|name={MHC_CHAIN_NAME}\n{mhc_seq.strip()}\n")
        fh.write(f">protein|name={PEPTIDE_CHAIN_NAME}\n{peptide_seq.strip()}\n")
    return fasta_path


def run_chai1_msa(
    peptide_seq: str,
    mhc_seq: str,
    output_dir: str | Path,
    *,
    num_trunk_samples: int = 5,
    num_diffn_samples: int = 5,
    num_diffn_timesteps: int | None = None,
    num_trunk_recycles: int | None = None,
    seed: int | None = None,
) -> Path:
    """
    Run Chai-1 on one pMHC complex with MSA ON (remote MMSeqs2/ColabFold
    server) and templates HARD-DISABLED.

    Returns the path to the Chai-1 structure/output directory
    (output_dir / "structures").
    """
    output_dir = Path(output_dir)
    fasta_path = build_fasta(peptide_seq, mhc_seq, output_dir)
    struct_dir = output_dir / "structures"

    # --- MSA-ON / TEMPLATE-OFF CONFIGURATION -- DO NOT CHANGE THE
    # TEMPLATE FLAGS WITHOUT RE-REVIEW ---
    # --use-msa-server: sequence-derived signal only, not a structure leak.
    # --no-use-templates-server, and --template-hits-path is NEVER set:
    # templates hand the model literal 3-D coordinates and are the actual
    # leak vector for this benchmark's post-2024-01-01 holdout.
    cmd = [
        CHAI_LAB_BIN,
        "fold",
        "--use-msa-server",
        "--no-use-templates-server",
        "--num-trunk-samples", str(num_trunk_samples),
        "--num-diffn-samples", str(num_diffn_samples),
        str(fasta_path),
        str(struct_dir),
    ]
    # -----------------------------------------------------------------

    if num_diffn_timesteps is not None:
        cmd += ["--num-diffn-timesteps", str(num_diffn_timesteps)]
    if num_trunk_recycles is not None:
        cmd += ["--num-trunk-recycles", str(num_trunk_recycles)]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    print(f"[run_chai1_msa] running: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)
    return struct_dir


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peptide-seq", required=True)
    p.add_argument("--mhc-seq", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--num-trunk-samples", type=int, default=5)
    p.add_argument("--num-diffn-samples", type=int, default=5)
    p.add_argument("--num-diffn-timesteps", type=int, default=None)
    p.add_argument("--num-trunk-recycles", type=int, default=None)
    p.add_argument("--seed", type=int, default=None)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_chai1_msa(
        peptide_seq=args.peptide_seq,
        mhc_seq=args.mhc_seq,
        output_dir=args.output_dir,
        num_trunk_samples=args.num_trunk_samples,
        num_diffn_samples=args.num_diffn_samples,
        num_diffn_timesteps=args.num_diffn_timesteps,
        num_trunk_recycles=args.num_trunk_recycles,
        seed=args.seed,
    )
