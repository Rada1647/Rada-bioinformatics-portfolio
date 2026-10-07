# Data sources and reconstruction

The release contains benchmark code, numerical results, and provenance records. It does not bundle simulation trajectories, topology files, extracted coordinates, or contact arrays. Reproduction downloads the original data from the providers and reconstructs the benchmark inputs locally. See [REPRODUCING.md](REPRODUCING.md) for commands.

The reference labels are contacts calculated from uncorrupted molecular dynamics simulations. They are not experimentally measured biological truth. All timing, construct identities, qualification decisions, and licenses below refer to the source versions inspected for this benchmark.

| Dataset | Benchmark role | Source version | Source-data license recorded at acquisition |
|---|---|---|---|
| T4 lysozyme | Primary test: one discovery run, one validation run, three held-out runs | [Zenodo 3989057](https://doi.org/10.5281/zenodo.3989057) | No explicit dataset license in the inspected metadata |
| Villin HP35 | Separate adaptive-sampling stress test | [Figshare 32541291, version 1](https://doi.org/10.6084/m9.figshare.32541291.v1) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| Ubiquitin | Qualification only; excluded from temporal denoising | [Zenodo 7792288](https://doi.org/10.5281/zenodo.7792288) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |

## T4 lysozyme

Credit: Felix Kümmerer, Simone Orioli, David Harding-Larsen, Falk Hoffmann, Yulian Gavrilov, Kaare Teilum, and Kresten Lindorff-Larsen. Dataset: *5x 1 µs all-atom MD trajectories; AMBER ff99SB*-ILDN & TIP4P/2005; T4 Lysozyme; Fitting side-chain NMR relaxation data using molecular simulations* (2020), [10.5281/zenodo.3989057](https://doi.org/10.5281/zenodo.3989057). Associated article: *Fitting Side-Chain NMR Relaxation Data Using Molecular Simulations*, JCTC 17(8), 5262–5275 (2021), [10.1021/acs.jctc.0c01338](https://doi.org/10.1021/acs.jctc.0c01338).

The deposit describes five independent 1 µs simulations. The protein-only XTC files contain 2,612 atoms and 1,000,001 frames at 1 ps intervals. Their matching TPR files describe the full 51,932-atom simulation systems. The benchmark extracts all 162 C-alpha atoms at exactly 200 ps intervals, giving 5,001 frames per run. It reconstructs the bonded protein across periodic boundaries before calculating distances.

The exact deposited construct is not assumed to be wild type. Relative to the inspected UniProt P00720 reference it has R12G, C54T, C97A, I137R, and lacks terminal NL. A reference difference alone does not establish that a substitution was engineered.

| Source trajectory | Role | Published and verified MD5 |
|---|---|---|
| `md1us4.xtc` | Discovery | `991c47f3ef4af3f0803a712264452dd6` |
| `md1us5.xtc` | Validation | `1a595d1027b671a3c467e70983151cc7` |
| `md1us6.xtc` | Test | `75807ff000533ee6f77f789dd2b4874f` |
| `md1us7.xtc` | Test | `1bc5ed165677a75d255528478bd59614` |
| `md1us8.xtc` | Test | `90c65520164769655aebd92d1760662d` |

The five XTCs total 50,659,155,984 bytes; the five matching TPRs bring the source download to 50,670,653,184 bytes, about 50.7 GB. Published sizes and topology checksums are retained in the source metadata. Complete-file MD5 verification is required after each download, including resumed downloads; local SHA-256 values are also recorded.

The source record marked the files open access but did not supply an explicit dataset license in the inspected metadata. This release therefore excludes both original T4 structural files and derived coordinate/contact arrays. It does not infer permissions from the separate ABSURDer software license. Download instructions do not assign a new license to the data.

## Villin HP35

Credit: Stefan Doerr, *HTMD tutorial data* (2026), Figshare version 1, [10.6084/m9.figshare.32541291.v1](https://doi.org/10.6084/m9.figshare.32541291.v1). Associated framework article: Stefan Doerr, Matthew J. Harvey, Frank Noé, and Gianni De Fabritiis, *HTMD: High-Throughput Molecular Dynamics for Molecular Discovery*, JCTC 12(4), 1845–1852 (2016), [10.1021/acs.jctc.6b00049](https://doi.org/10.1021/acs.jctc.6b00049). Source conditions also come from the [official villin-folding tutorial](https://software.acellera.com/htmd/tutorials/analysis/villin-folding.html).

- File: `protein_folding_datasets.zip`, [direct source download](https://ndownloader.figshare.com/files/65180772).
- Size: 2,559,410,654 bytes.
- Published and verified MD5: `484eafab71b27e3d66f48d4820d08f47`.
- Local SHA-256: `4f512459ac7bfc42e8bcd6b247648b9ccdfe60dde7d8ace1dbeb765293618fff`.
- Source metadata: [Figshare API record](https://api.figshare.com/v2/articles/32541291).

The archive contains 2,137 trajectories in three groups of 708, 710, and 719. The source tutorial describes 360 K simulations with 100 ps frame spacing. Extraction retains every second frame, preserves trajectory boundaries and recorded source times, and selects all 35 C-alpha atoms, including the noncanonical residues. The topology contains NLE65, HSP68, and NLE70; this is an engineered HP35 construct, not wild-type villin. PDB [2F4K](https://www.rcsb.org/structure/2F4K) corroborates the residue pattern, but does not establish the simulation's exact starting structure or pH.

The frozen protocol assigns groups 1, 2, and 3 to discovery, validation, and test. Within each group, ascending SHA-256 hashes of canonical source trajectory paths select the first 24, 24, and 48 eligible trajectories. All 2,137 trajectories qualified; no trajectory was excluded using denoising outcomes. One of the 48 test trajectories has 71 sampled frames and the other 47 have 250.

Adaptive ancestry makes these trajectories dependent. The 48 test paths include six ancestral roots within one held-out group; they are not 48 independent biological replicates. Independence between the three groups is not certified by the available source metadata. Villin is therefore reported separately and cannot change the primary T4 decision.

The source dataset is explicitly CC BY 4.0. This benchmark changes the representation by subsampling time points, extracting C-alpha coordinates, selecting discovery-variable contact features, introducing specified artificial corruptions, and computing aggregate scores. The original coordinates and derived arrays are omitted from this compact release; metadata and attribution remain. No endorsement by the source authors is implied.

## Ubiquitin

Credit: Adekunle Aina, Shawn C. C. Hsueh, and Steven S. Plotkin, *Molecular dynamics-generated ensemble dataset of ubiquitin; for PROTHON: A Local Order Parameter-Based Method for Efficient Comparison of Protein Ensembles* (2023), [10.5281/zenodo.7792288](https://doi.org/10.5281/zenodo.7792288). Associated article: *PROTHON: A Local Order Parameter-Based Method for Efficient Comparison of Protein Ensembles*, JCIM 63(11), 3453–3461 (2023), [10.1021/acs.jcim.3c00145](https://doi.org/10.1021/acs.jcim.3c00145).

- File: `ubiquitin-md-generated-ensemble.zip`, 228,792,237 bytes.
- Published and verified MD5: `3372f857beb371ca5ebc8d57877c5c79`.
- Local SHA-256: `d21f88be709745f858a78c11dbfb7ff6bdb014e6cb0aaac2b489071f540dd0d8`.
- License in the saved Zenodo record: CC BY 4.0.

The archive contains a 76-residue topology and six 5,000-frame DCDs. Q99 is the native reference ensemble; Q75–Q95 are biased partially folded ensembles. These are not six independent native simulations.

The article describes five native runs, but the inspected archive and author repositories did not establish the mapping from archived Q99 frame indices to original run boundaries. The DCD-header interval interpreted by MDAnalysis also disagreed with the article's sampling interval. Splitting 5,000 frames into five consecutive blocks of 1,000 would be an unsupported assumption.

The predefined temporal-provenance gate therefore failed. No ubiquitin contact-denoising scores were computed. This is a limitation of this benchmark's available provenance, not a finding that the published structural ensemble is invalid. The inspection was bounded to the documented sources; inaccessible supplementary material was not claimed to have been inspected, and the source authors were not contacted. Reliable assembly metadata or original separate trajectories would be needed to qualify it.

## Provenance and audit boundaries

Historical metadata and extraction manifests are evidence about the original run, not generated inputs for a new run. Keeping them in a provenance directory prevents completion markers from making acquisition skip data that have not been reconstructed locally. Newly generated manifests can contain different timestamps or paths even when their scientific arrays match.

The protocol was privately frozen before new-protein outcomes. It was not publicly preregistered. Its authoritative v1.1 JSON has SHA-256 `1e25ff2c9890809f5432675c259f8ec705d2ce9bcca23969111f771fd520c84c`. Public-package hashes describe the released bytes; original source and computation hashes describe their respective original artifacts. A redacted or adapted release file must not be represented as byte-identical to its original.

Five independent AI review tasks examined mathematics, data, statistics, citations, and reproduction. They are technical audits, not human peer review. Full numerical reruns used a fresh benchmark environment on the same host, not an independent machine. For T4, the final raw trajectory received an additional independent source-level audit; the other four raw files were not independently downloaded again after extraction. These boundaries are preserved in the review records.
