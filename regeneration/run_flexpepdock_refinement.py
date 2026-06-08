"""Run FlexPepDock refinement on all threaded inputs from the regeneration.

Locates the {peptide}_input.pdb files produced by the threading step,
then runs FlexPepDocking.linuxgccrelease with -nstruct 25 to produce a
25-decoy ensemble + score.sc per peptide–allele pair. Output is
organized into the same layout the validation notebook expects:

    regeneration/pdb/{allele_dir}/{peptide}/
        {peptide}_input_0001.pdb
        ...
        {peptide}_input_0025.pdb
        score.sc

Parallelizes across peptides with multiprocessing.Pool; each FlexPepDock
run is itself single-threaded.

Usage:
    python run_flexpepdock_refinement.py
"""

from __future__ import annotations

import multiprocessing as mp
import shutil
import subprocess
import time
from pathlib import Path

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

REGEN_BASE       = Path('/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration')
THREADING_OUTPUT = REGEN_BASE / 'fastas' / 'output'
PDB_OUTPUT       = REGEN_BASE / 'pdb'
DOCKING_LOGS     = REGEN_BASE / 'docking_logs'

ROSETTA_BIN = '/sb/meilerapps/rosetta/rosetta-3.15/main/source/bin/FlexPepDocking.linuxgccrelease'
ROSETTA_DB  = '/sb/meilerapps/rosetta/rosetta-3.15/main/database'

NSTRUCT     = 25
N_PARALLEL  = 28


# -----------------------------------------------------------------------------
# Discovery
# -----------------------------------------------------------------------------

def find_threaded_inputs() -> list[dict]:
    """Locate every {peptide}_input.pdb under the threading output tree.

    Directory layout produced by IEDBTestPipeline_ACCRE.py:
        fastas/output/{ALLELE_DIR}_validation_batch{N}/{PEPTIDE}/{PEPTIDE}_input.pdb
    """
    pairs = []
    for input_pdb in THREADING_OUTPUT.rglob('*_input.pdb'):
        peptide = input_pdb.parent.name
        batch_dir = input_pdb.parent.parent.name  # e.g. 'A0101_validation_batch1'
        allele_dir = batch_dir.split('_validation_')[0]
        pairs.append({
            'allele_dir': allele_dir,
            'peptide':    peptide,
            'input_pdb':  input_pdb,
        })
    return pairs


# -----------------------------------------------------------------------------
# Per-pair worker
# -----------------------------------------------------------------------------

def run_one(pair: dict) -> dict:
    """Run FlexPepDock refinement for one peptide–allele pair."""
    allele_dir = pair['allele_dir']
    peptide    = pair['peptide']
    input_pdb  = pair['input_pdb']

    out_dir = PDB_OUTPUT / allele_dir / peptide
    out_dir.mkdir(parents=True, exist_ok=True)

    # Copy the threaded input into the output dir so Rosetta-named decoys
    # come out as {peptide}_input_NNNN.pdb (matching the released dataset).
    local_input = out_dir / input_pdb.name
    if not local_input.exists():
        shutil.copy2(input_pdb, local_input)

    cmd = [
        ROSETTA_BIN,
        '-in:file:s',          str(local_input),
        '-database',           ROSETTA_DB,
        '-pep_refine',
        '-nstruct',            str(NSTRUCT),
        '-ex1',
        '-ex2aro',
        '-out:file:scorefile', 'score.sc',
        '-overwrite',
    ]

    DOCKING_LOGS.mkdir(parents=True, exist_ok=True)
    log_path = DOCKING_LOGS / f'{allele_dir}_{peptide}.log'

    t0 = time.time()
    with open(log_path, 'w') as logf:
        result = subprocess.run(cmd, cwd=str(out_dir), stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0

    # Quick sanity check on output
    decoys   = sorted(out_dir.glob(f'{peptide}_input_*.pdb'))
    score_sc = out_dir / 'score.sc'
    n_decoys = len(decoys)
    has_score = score_sc.is_file()

    return {
        'peptide':      peptide,
        'allele_dir':   allele_dir,
        'returncode':   result.returncode,
        'n_decoys':     n_decoys,
        'has_score_sc': has_score,
        'elapsed_min':  round(elapsed / 60, 1),
        'log_path':     str(log_path),
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    pairs = find_threaded_inputs()
    #pairs = find_threaded_inputs()[:1] # To test first allele
    print(f'Found {len(pairs)} threaded inputs to refine')
    print(f'Output base:  {PDB_OUTPUT}')
    print(f'Docking logs: {DOCKING_LOGS}')
    print(f'Parallelism:  {N_PARALLEL} workers')
    print(f'Per pair:     {NSTRUCT} decoys × ~3 min ≈ 75 min single-threaded')
    print()
    print(f'Expected wall time: ~{(len(pairs) * 75) // N_PARALLEL} min '
          f'(plus startup overhead)')
    print()

    t_start = time.time()
    with mp.Pool(N_PARALLEL) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(run_one, pairs), 1):
            results.append(result)
            status = 'OK ' if (result['returncode'] == 0 and result['n_decoys'] == NSTRUCT) else 'FAIL'
            print(f'[{i:3d}/{len(pairs)}] {status} '
                  f'{result["allele_dir"]}/{result["peptide"]}  '
                  f'({result["n_decoys"]}/{NSTRUCT} decoys, '
                  f'{result["elapsed_min"]} min)')

    elapsed_total = time.time() - t_start
    print()
    print('=' * 70)
    print(f'Total time: {elapsed_total/60:.1f} min')
    success = sum(1 for r in results if r['returncode'] == 0 and r['n_decoys'] == NSTRUCT)
    print(f'Successful: {success}/{len(pairs)}')
    fail = [r for r in results if not (r['returncode'] == 0 and r['n_decoys'] == NSTRUCT)]
    if fail:
        print(f'Failed:     {len(fail)}')
        for r in fail:
            print(f'  {r["allele_dir"]}/{r["peptide"]}: '
                  f'rc={r["returncode"]}, n_decoys={r["n_decoys"]}, '
                  f'log={r["log_path"]}')


if __name__ == '__main__':
    main()
