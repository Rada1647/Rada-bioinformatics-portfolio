# Reproducing the negative result

This repository contains the frozen protocol, benchmark code, numerical results,
source metadata and review records. It does **not** contain protein coordinates,
contact arrays, topologies or source trajectory archives. Those inputs must be
reconstructed locally before a protein benchmark rerun. See
[DATA_SOURCES.md](DATA_SOURCES.md) for the pinned datasets and rights information.

Run all commands from the repository root. The tested environment was Linux,
Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. Source extraction additionally used
MDAnalysis 2.9.0. The T4 downloader uses POSIX `os.pwrite`; use Linux or a suitable
Linux environment for the documented full reconstruction.

## 1. Fast offline checks

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-benchmark.txt
.venv/bin/python scripts/check_release.py
```

This checks release file hashes, the resolvent against a dense inverse, the cubic
against explicit expected-median enumeration, and a synthetic end-to-end
benchmark including output-overwrite protection. It also rebuilds `summary.json`
and `metrics.csv` from archived numerical results in a temporary directory.
No protein data are downloaded and no protein benchmark is rerun by these checks.

To inspect only report aggregation, Python's standard library is sufficient:

```bash
python3 scripts/rebuild_summary.py
python3 scripts/rebuild_summary.py --out local_outputs/summary
```

The second command retains regenerated files in a **new** directory. Scientific
summary content must match exactly after excluding only `generated_utc`;
`metrics.csv` must match byte for byte. This verifies aggregation of stored
numbers, not independent correctness of their underlying benchmark scores.

## 2. Reconstruct T4 lysozyme from original deposited trajectories

```bash
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python data/t4/acquire_extract.py --workers 16
.venv/bin/python prepare_extension.py --dataset t4
```

The acquisition script reads the saved Zenodo record for
[10.5281/zenodo.3989057](https://doi.org/10.5281/zenodo.3989057). It checks the
complete published size and MD5 for each downloaded file, records SHA256, and
extracts 162 C-alpha coordinates at exactly 200 ps. Bonded backbone reconstruction
handles periodic boundaries. All five runs are required:

| Source run | Frozen use |
|---|---|
| `md1us4` | Discovery of variable contacts |
| `md1us5` | Hyperparameter validation |
| `md1us6`, `md1us7`, `md1us8` | Held-out primary test runs |

T4 acquisition transfers about 51 GB. Each raw trajectory is about 10 GB. The
script processes runs sequentially and removes newly downloaded raw XTCs for
runs 4–7 only after checksum and extraction succeed; it retains `md1us8.xtc`.
`--keep-raw` retains all newly downloaded trajectories. Files present before the
invocation are never removed by this policy. Allow ample disk space for the raw
file in progress, extracted arrays and any retained files; the two full datasets
together require roughly 54 GB of transfer, and a 25 GB free-space budget is a
reasonable starting allowance with the default cleanup behavior. Retaining every
raw T4 trajectory requires more than 51 GB of disk space.

The public-release extractor adds one operational guard to the original: a
completion record causes a skip only if its coordinate output exists and matches
the recorded size and SHA256. Scientific extraction operations are unchanged.
Historical completion records are in `provenance/t4/`, not the live `data/t4/`
output directory. Do not copy those records into live output directories as a
substitute for extraction. A stale or mismatched completion record fails clearly;
preserve the problematic files for inspection and retry from a fresh checkout.

Preparation enforces provenance and geometry checks, creates
`data/prepared/t4/`, and refuses to overwrite an existing prepared directory.
Contact feature selection uses discovery data only. The checksum-verified raw
source is the reference; this task does not rerun the underlying molecular
dynamics simulations.

## 3. Reconstruct Villin under its separate frozen protocol

```bash
.venv/bin/python scripts/download_source.py villin
.venv/bin/python data/villin/qualify_extract.py
.venv/bin/python prepare_extension.py --dataset villin
```

The downloader checks the exact byte count, published MD5 and recorded SHA256 of
the 2.56 GB archive from
[10.6084/m9.figshare.32541291.v1](https://doi.org/10.6084/m9.figshare.32541291.v1).
The extractor, copied unchanged from the experiment, checks the saved Figshare
metadata and preserves all 2,137 trajectory boundaries. It samples every second
100 ps source frame, without interpolating or joining paths.

Preparation selects 24 discovery, 24 validation and 48 test trajectories from
archive groups 1, 2 and 3 respectively, ordered by SHA256 of the canonical source
path after the predefined metadata eligibility checks. The complete manifest
order determines each trajectory's random-seed index. Do not replace it with
filesystem enumeration or select trajectories by their observed contacts.

Villin is an adaptive-sampling **stress test**, not an additional independent
confirmatory protein experiment. Archive-group independence is unverified;
the 48 held-out paths are not 48 independent biological replicates.

## 4. Rerun all frozen benchmark scenarios

After reconstructing the relevant dataset, use new output directories:

```bash
.venv/bin/python run_scenarios.py --dataset t4 --workers 4 --out local_outputs/t4
.venv/bin/python run_scenarios.py --dataset villin --workers 4 --out local_outputs/villin
.venv/bin/python compare_reproduction.py results/t4 local_outputs/t4 --allow-regenerated-manifest --out local_outputs/t4_comparison.json
.venv/bin/python compare_reproduction.py results/villin local_outputs/villin --allow-regenerated-manifest --out local_outputs/villin_comparison.json
```

Use `--workers 1` if memory is constrained. This changes scheduling, not the fixed
random seeds. The scenario runner executes all eight corruption/cutoff scenarios;
the primary scenario also writes its clean-input control. Each dataset should
produce eight validation JSONs, eight test JSONs and `primary_clean.json`.

Do not use `results/t4` or `results/villin` as new output targets: those are the
archived originals. Completed scientific outputs are protected against overwrite.
If an invocation fails after producing partial outputs, retain it for diagnosis
and rerun into a new directory.

`--allow-regenerated-manifest` permits the preparation manifest hash to differ
because extraction provenance and preparation timestamps are new. It does **not**
waive prepared array hashes, protocol or algorithm hashes, feature identities,
parameters, scores or decisions. Other excluded metadata fields are completion
times, runtimes and the validation file hash, which itself changes with those
timestamps. The comparison script reports its exact exclusions.

Exact equality is the strict reproduction target; floating-point libraries or
different extraction versions may introduce differences that require inspection.
The historical fresh-environment check reran the benchmark on the same host with
the same prepared arrays. It did not establish an independent-machine, fresh
raw-data reconstruction. The independent reviews record their separate scopes.
Threshold classification can be sensitive to roundoff extremely close to 0.5;
the audited cases did not change the primary failure decision. Do not automatically
discard discrepancies as harmless roundoff.

## 5. Reinspect the excluded ubiquitin archive

```bash
.venv/bin/python scripts/qualify_ubiquitin.py
```

This downloads or verifies the pinned 229 MB archive from
[10.5281/zenodo.7792288](https://doi.org/10.5281/zenodo.7792288), checks ZIP members,
extracts locally, and writes `data/ubiquitin/archive_metadata_reproduced.json`.
It reports topology counts, six 5,000-frame DCD files and their header timing.
For a repeat inspection, preserve or move the existing `data/ubiquitin/archive/`
directory first; the helper deliberately refuses to reuse it.

Compare these findings with `provenance/ubiquitin/archive_metadata.json`,
`qualification.json` and `source_provenance.json`. The helper reproduces archive
inspection only. It does not repeat the historical paper/repository search or
establish that no mapping exists anywhere. The recorded search did not establish
original run-to-frame mapping, and MDAnalysis's DCD interval of 0.04888821 ps did
not match the paper's reported 40 ps sampling. No five-by-1,000-frame partition was
assumed. No ubiquitin contact-denoising benchmark was performed or is enabled by
these reproduction commands.

## What would count as the same result?

The primary T4 decision must remain `FAIL`: the cubic does not meet the frozen
combined practical-success rule against the validation-selected strongest
prespecified nonoracle control. The matched seven-setting power control also wins.
Villin's separately reported stress criterion fails. Reproducing these conclusions
does not disprove the source theorem; the theorem does not promise superiority in
this Brier-error denoising task.

The exact archived protocol is `protocol/protocol_v1.1.json`, SHA256
`1e25ff2c9890809f5432675c259f8ec705d2ce9bcca23969111f771fd520c84c`.
It was frozen privately before the new-protein outcomes, not publicly preregistered.
Benchmark algorithms and the protocol remain unchanged in this release. Added
release helpers arrange downloads and checks; they do not tune methods or alter
scientific outputs.
