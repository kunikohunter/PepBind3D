"""
Filesystem roots, in one place and overridable by environment variable.

Every script in this repository derives its paths from here rather than
hardcoding them, so running the code somewhere else means setting two
environment variables instead of editing 26 files:

    export PEPBIND3D_DATA=/your/path/to/IEDB_data_clean
    export PEPBIND3D_MHC_DB=/your/path/to/MHC_database
    export PEPBIND3D_CLUSTER=/your/cluster/scratch

The defaults are the authors' cluster layout, so the code reproduces the
published analyses unchanged when run in place.

Deliberately pure stdlib and importable without side effects: the silent-file
converter runs under PyRosetta fanned out over 112k pairs, and `utils/__init__`
eagerly imports biopython and matplotlib, so this cannot live there.
"""
import os
from pathlib import Path

# Curated measurements, structural ensembles and all analysis outputs.
#   metadata.csv, structures/, release_v2_final/, IEDB_validation/, pdb/
DATA_ROOT = Path(os.environ.get(
    "PEPBIND3D_DATA",
    "<HOME>/main_project/data/IEDB_data_clean"))

# Local MHC template database: threading templates, reference crystals,
# database.info, default_receptor/, and the raw IEDB bulk export under build/.
MHC_DB_ROOT = Path(os.environ.get(
    "PEPBIND3D_MHC_DB",
    "<HOME>/Data/MHC_database"))

# Cluster scratch root for run trees that live on the allocation. Nothing in
# the release pipeline or the analyses needs it.
CLUSTER_ROOT = Path(os.environ.get(
    "PEPBIND3D_CLUSTER",
    "<CLUSTER>"))

# Rosetta installation root. Only needed to generate structures or to run
# extract_pdbs; the analysis scripts do not use it.
ROSETTA_ROOT = Path(os.environ.get(
    "PEPBIND3D_ROSETTA",
    "/dors/meilerlab/apps/rosetta/rosetta-3.14/main"))

# Frequently used files, named so a reader does not have to reconstruct them.
TEMPLATE_DIR = MHC_DB_ROOT / "templates"
MHC_DB_INFO = MHC_DB_ROOT / "database.info"
RAW_IEDB = MHC_DB_ROOT / "build" / "mhc_ligand_full.csv"
RELEASE_DIR = DATA_ROOT / "release_v2_final"
RELEASE_METADATA = RELEASE_DIR / "metadata.csv"
RELEASE_STRUCTURES = RELEASE_DIR / "structures"
VALIDATION_DIR = DATA_ROOT / "IEDB_validation"


def describe():
    """Print the resolved roots and whether they exist. Useful as a first check
    when running this code on a new machine."""
    for name, p in [("DATA_ROOT", DATA_ROOT), ("MHC_DB_ROOT", MHC_DB_ROOT),
                    ("CLUSTER_ROOT", CLUSTER_ROOT), ("ROSETTA_ROOT", ROSETTA_ROOT)]:
        env = {"DATA_ROOT": "PEPBIND3D_DATA", "MHC_DB_ROOT": "PEPBIND3D_MHC_DB",
               "CLUSTER_ROOT": "PEPBIND3D_CLUSTER",
               "ROSETTA_ROOT": "PEPBIND3D_ROSETTA"}[name]
        src = "env" if env in os.environ else "default"
        print(f"  {name:13s} {str(p):55s} [{src}] "
              f"{'ok' if p.exists() else 'MISSING'}")


if __name__ == "__main__":
    describe()
