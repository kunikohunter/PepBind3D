import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

PDB_DIR = Path("/home/huntek1/main_project/data/IEDB_data_clean/pdb")
HF_STRUCT_DIR = Path("/home/huntek1/main_project/data/IEDB_data_clean/huggingface/structures")

# Columns in score.sc's "description" field are the decoy tag, not a score --
# everything else gets carried into the silent file as-is via
# setPoseExtraScore, so the FlexPepDock interface terms (I_sc, I_bsa, I_hb,
# I_pack, I_unsat, pep_sc, rmsBB, ...) survive instead of being replaced by a
# generic score_jd2 rescore.
TAG_COLUMN = "description"


def parse_score_sc(score_sc_path):
    """Parse a Rosetta score.sc into {tag: {column: float}}. Returns {} if
    the file is missing or has no SCORE: data rows (caller decides how to
    treat that)."""
    rows = {}
    header = None
    with open(score_sc_path) as f:
        for line in f:
            if not line.startswith("SCORE:"):
                continue
            parts = line.split()[1:]  # drop the leading "SCORE:" token
            if header is None:
                header = parts  # first SCORE: line is the column header
                continue
            if len(parts) != len(header):
                continue  # malformed row; skip rather than misalign columns
            row = dict(zip(header, parts))
            tag = row.pop(TAG_COLUMN)
            scores = {}
            for k, v in row.items():
                try:
                    scores[k] = float(v)
                except ValueError:
                    pass  # non-numeric column; skip
            rows[tag] = scores
    return rows


_pyrosetta_ready = False


def _worker_init():
    """Runs once per ProcessPoolExecutor worker (not once per peptide) --
    PyRosetta init is expensive and must happen in each worker process."""
    global _pyrosetta_ready
    import pyrosetta
    pyrosetta.init("-mute all")
    _pyrosetta_ready = True


def convert_peptide(args):
    allele, peptide, pep_dir = args
    pep_dir = Path(pep_dir)
    out_dir = HF_STRUCT_DIR / allele
    out_dir.mkdir(parents=True, exist_ok=True)
    silent_out = out_dir / f"{peptide}.silent"

    if silent_out.exists() and silent_out.stat().st_size > 0:
        return allele, peptide, "skipped"

    score_sc = pep_dir / "score.sc"
    if not score_sc.exists():
        return allele, peptide, "no_score_sc"
    scores_by_tag = parse_score_sc(score_sc)
    if not scores_by_tag:
        return allele, peptide, "empty_score_sc"

    # Only numbered decoy PDBs -- explicitly excludes any pre-dock/relaxed
    # starting model file (e.g. "{peptide}_input.pdb", no number suffix) that
    # may sit alongside the numbered decoys in some peptide directories.
    pdbs = sorted(pep_dir.glob(f"{peptide}_input_[0-9][0-9][0-9][0-9].pdb"))
    if not pdbs:
        return allele, peptide, "no_pdbs"

    import pyrosetta
    from pyrosetta.rosetta.core.io.silent import SilentFileData, SilentFileOptions, BinarySilentStruct
    from pyrosetta.rosetta.core.pose import setPoseExtraScore

    opts = SilentFileOptions()
    sfd = SilentFileData(opts)

    n_written = 0
    tags_missing_scores = []
    for pdb_path in pdbs:
        # Two DIFFERENT identifiers, deliberately kept separate:
        #  - score_key matches score.sc's "description" column, which has no
        #    trailing suffix (e.g. "AADFPGIAR_input_0001").
        #  - tag is what goes into the silent file. v1 tags carry a trailing
        #    "_0001", an artifact of the old per-file score_jd2/JD2 conversion
        #    (e.g. "AADFPGIAR_input_0001_0001"). It is preserved so regenerated
        #    files stay tag-compatible with the published release; anything
        #    keyed on v1 tags would silently stop matching otherwise.
        score_key = pdb_path.stem
        tag = f"{score_key}_0001"
        row = scores_by_tag.get(score_key)
        if row is None:
            tags_missing_scores.append(score_key)
            continue
        pose = pyrosetta.pose_from_pdb(str(pdb_path))
        for col, val in row.items():
            setPoseExtraScore(pose, col, val)
        ss = BinarySilentStruct(opts)
        ss.fill_struct(pose, tag)
        sfd.add_structure(ss)
        n_written += 1

    if n_written == 0:
        return allele, peptide, "no_decoys_written"

    # Data-integrity check: every PDB found on disk must have made it into
    # the silent file. The old script's failure mode (25 PDBs in, 24 decoys
    # out, silently) is exactly what this catches -- report it as an ERROR
    # instead of writing a silently-incomplete file.
    if n_written != len(pdbs):
        return allele, peptide, (
            f"ERROR: found {len(pdbs)} PDBs but only wrote {n_written} decoys "
            f"(missing score.sc rows for: {tags_missing_scores})"
        )

    sfd.write_all(str(silent_out))
    return allele, peptide, "ok"


def convert_peptide_safe(args):
    """Wraps convert_peptide so a worker-process crash on one peptide
    (e.g. a malformed PDB) doesn't kill the whole batch silently."""
    try:
        return convert_peptide(args)
    except Exception as e:
        allele, peptide, _ = args
        return allele, peptide, f"ERROR: {type(e).__name__}: {e}"


if __name__ == "__main__":
    tasks = []
    for allele_dir in sorted(PDB_DIR.iterdir()):
        if not allele_dir.is_dir():
            continue
        for pep_dir in sorted(allele_dir.iterdir()):
            if not pep_dir.is_dir():
                continue
            tasks.append((allele_dir.name, pep_dir.name, str(pep_dir)))

    print(f"Total peptide dirs to convert: {len(tasks)}")

    errors = []
    # initializer= runs pyrosetta.init() once per worker process, not once
    # per peptide -- PyRosetta init is too slow to redo for every task.
    with ProcessPoolExecutor(max_workers=16, initializer=_worker_init) as executor:
        futures = {executor.submit(convert_peptide_safe, t): t for t in tasks}
        with tqdm(total=len(tasks)) as pbar:
            for future in as_completed(futures):
                allele, peptide, status = future.result()
                if "ERROR" in str(status):
                    errors.append((allele, peptide, status))
                pbar.update(1)

    print(f"\nDone. {len(errors)} errors.")
    if errors:
        for e in errors:
            print(e)
