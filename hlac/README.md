# HLA-C generation (v2)

Generates FlexPepDock ensembles for 1,861 HLA-C peptide-allele pairs across
11 alleles, extending the released HLA-A/B dataset for v2.

## Why this is separate from `regeneration/`

The two directories look similar and do opposite things. Keep them apart.

|                  | `regeneration/`                        | `hlac/`                          |
|------------------|----------------------------------------|----------------------------------|
| Scope            | 52 already-released validation pairs    | 1,861 new pairs                  |
| Self-templates   | **excluded** (`--ignore_epitope_match`) | **allowed** (`omit=[]`)          |
| Rosetta          | `rosetta-3.15` tree                     | `2024.09+release.06b3cf8`        |
| Execution        | serial loop, 28 threads, one host       | SLURM arrays with `%` throttling |

`regeneration/` exists specifically because self-templating was a defect for
the validation subset: those pairs were threaded onto their own native
crystals, so the structural validation measured refinement perturbation
rather than modeling accuracy. Excluding self-templates is the entire point
of that directory.

HLA-C is the opposite case. Self-matches are allowed, because that is how the
released HLA-A/B structures were built, and v2 must stay internally
consistent with them. Running HLA-C with self-exclusion would silently break
that consistency and would not be visible in the output.

The Rosetta versions also differ. All three current arms (released dataset,
HLA-C, and the benchmark) stamp `REMARK 220 VERSION 2024.09+release.06b3cf8`;
`regeneration/` targets an older tree. Merging the two would force one of
them onto the wrong build.

## Workflow

```
make_threading_inputs.py     # 1,861 targets -> 81 per-allele batch FASTAs
jobs/thread_array.sbatch     # array 1-81, threads peptides onto receptors
jobs/stage.sbatch            # runs stage_and_verify.py, gates on exact count
jobs/dock_array.sbatch       # array 1-1861%100, FlexPepDock nstruct 25, 1G
```

`stage_and_verify.py` fails unless it finds exactly 1,861 threaded inputs.
This gate is load-bearing: `IEDBTestPipeline_ACCRE.py` uses `os.mkdir` rather
than `os.makedirs(exist_ok=True)`, and on a rerun over an existing directory
it raises `[Errno 17]`, swallows it, and exits 0. SLURM then reports the task
COMPLETED with structures missing. Never trust the exit code here; check the
count.

Resource requests match the production `flexpep_array_job.sh`: 1 GB and one
CPU per task. Docking measured ~75 min per pair for 25 decoys.

## Verification

Every threaded structure was checked against the receptor for its own allele:
chain A compared byte-for-byte with `default_receptor/{allele}.pdb`, chain B
against the peptide its directory is named for. All 1,861 passed, with
per-allele counts reproducing `hlac_threading_targets.csv` exactly.
