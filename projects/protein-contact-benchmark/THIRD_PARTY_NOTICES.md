# Third-party sources and notices

Any license applied to this repository's original code or writing does not relicense external datasets, papers, software dependencies, or upstream source records. This release contains no third-party article PDFs, cloned author repositories, molecular trajectory/topology files, or extracted coordinate/contact arrays.

## Source datasets

**T4 lysozyme.** Kümmerer, Orioli, Harding-Larsen, Hoffmann, Gavrilov, Teilum, and Lindorff-Larsen (2020), [Zenodo dataset 10.5281/zenodo.3989057](https://doi.org/10.5281/zenodo.3989057). No explicit dataset license was supplied in the metadata inspected for this benchmark. The repository includes citations, acquisition code, source metadata, checksums, and numerical analysis records; it excludes original and derived structural arrays. Open access was not treated as an explicit redistribution license. The associated [JCTC article](https://doi.org/10.1021/acs.jctc.0c01338) and separate ABSURDer software are not bundled or relicensed here.

**Villin HP35.** Stefan Doerr (2026), *HTMD tutorial data*, Figshare, version 1, [10.6084/m9.figshare.32541291.v1](https://doi.org/10.6084/m9.figshare.32541291.v1), licensed by its provider under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). The analysis subsamples trajectories, extracts C-alpha coordinates, builds contact features, adds artificial corruption, and reports derived scores. Source archives and coordinate/contact arrays are not bundled. Source metadata and provenance summaries retain this attribution. Also credit Doerr, Harvey, Noé, and De Fabritiis, [HTMD (2016)](https://doi.org/10.1021/acs.jctc.6b00049), and the [Acellera tutorial](https://software.acellera.com/htmd/tutorials/analysis/villin-folding.html).

**Ubiquitin.** Adekunle Aina, Shawn C. C. Hsueh, and Steven S. Plotkin (2023), *Molecular dynamics-generated ensemble dataset of ubiquitin; for PROTHON: A Local Order Parameter-Based Method for Efficient Comparison of Protein Ensembles*, [10.5281/zenodo.7792288](https://doi.org/10.5281/zenodo.7792288), licensed by its provider under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). Only qualification metadata and findings are reported; no denoising scores were computed. The associated [PROTHON article](https://doi.org/10.1021/acs.jcim.3c00145), original archive, and source repositories are not bundled.

The dataset license statements above describe the saved acquisition metadata. They do not assert ownership by this benchmark's author or endorsement by the data providers. Full dataset identities, transformations, and checksums are in [DATA_SOURCES.md](DATA_SOURCES.md).

## Mathematical and scientific attribution

The benchmark studies a construction discussed in OpenAI's [*Metric Markov Cotype Two of ℓ1*, October 5, 2026](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Metric-Markov-Cotype-Two-of-l1-October-5-2026/l1-markov-cotype.pdf), pinned to repository commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`. The binary cubic `3h² − 2h³` and nonlinear cut smoothing have earlier attribution to Qingjin Cheng, Yue Wang, and Bo Xiang, [*Sharp Metric Cotype Inequalities for L1 via Nonlinear Cut Smoothing*, arXiv:2609.08749v2](https://arxiv.org/abs/2609.08749v2). The cubic itself is not presented as newly invented by this benchmark or by the later proof.

The mathematical source motivates the tested construction. Its inequality does not guarantee superior Brier denoising scores, recovery of protein kinetics, or physically valid protein structures. Failure of the engineering benchmark is not a disproof of that theorem. Source papers are cited, not redistributed. The complete study bibliography and claim scopes are retained in `references.json` and the accompanying report.

## Software dependencies

Python, NumPy, SciPy, MDAnalysis, Matplotlib, and ReportLab are installed separately under their own licenses. Their source code and distribution packages are not vendored in this release. Version pins describe the recorded computational environment; they do not override dependency licenses. Scientific software citations, including both MDAnalysis references, appear in `references.json`.

No association with, approval by, or endorsement from OpenAI, the cited researchers, the data providers, or the dependency maintainers is claimed.
