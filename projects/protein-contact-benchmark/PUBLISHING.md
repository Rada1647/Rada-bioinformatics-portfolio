# Location in the portfolio

This benchmark is maintained at:

https://github.com/Rada1647/Rada-bioinformatics-portfolio/tree/main/projects/protein-contact-benchmark

Run reproduction commands from this project directory, following REPRODUCING.md.
The active automated check is `.github/workflows/benchmark-checks.yml` at the
portfolio root; it installs the pinned benchmark requirements and invokes
`scripts/check_release.py` in this project. Source datasets are downloaded
separately and must not be added to the repository.

The scientific code, frozen protocol, original numerical results and audit
provenance remain the archived study. Later experiments should use a separate
protocol and output directory. Project documentation can evolve while preserving
the identity of those scientific files and recording release-manifest updates.

No project DOI has been minted. The citation file points to the actual repository.
