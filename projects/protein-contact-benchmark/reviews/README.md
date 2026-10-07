# Technical review records

Five independent AI review roles examined the completed benchmark:

| Record | Scope |
|---|---|
| [01_math_final.json](01_math_final.json) | Resolvent/cubic mathematics, implementation and independent numerical checks |
| [02_data_final.json](02_data_final.json) | Source qualification, trajectory extraction, checksums and contact construction |
| [03_statistics_final.json](03_statistics_final.json) | Frozen rules, selections, aggregation, uncertainty and reported numbers |
| [04_citations_final.json](04_citations_final.json) | References, attribution and scope of scientific claims |
| [05_reproduction_final.json](05_reproduction_final.json) | Fresh-environment reruns and alternative numerical implementations |

These are **AI technical audits, not human peer review**. Each file records its
own scope and limitations. See [REVIEW_RESPONSE.txt](../REVIEW_RESPONSE.txt) for
findings and resolutions. Supporting JSON and historical audit scripts accompany
the final records; some require source arrays, the original private package, or
paths from the original environment and are not turnkey public entry points.
The maintained public entry point is [REPRODUCING.md](../REPRODUCING.md).

The statistics review is an explicitly redacted copy: private storage identifiers
and local storage metadata were removed, with original/release hashes and exact
redacted JSON paths in [RELEASE_ORIGINS.json](../provenance/RELEASE_ORIGINS.json).
All other included historical review files are unchanged. Historical mentions
of unavailable private receipts are not evidence that those receipts are included.

The [T4](fresh_t4_comparison.json) and [villin](fresh_villin_comparison.json)
comparison reports document exact reruns across 17 output files per dataset.
The reruns used newly installed benchmark dependencies on the same host. They
are not independent-computer or independent-laboratory replication. Absolute
scratch paths in historical records describe that past run and are not paths
that a reproducer must create.

Alternative spectral implementations found a small number of scores close to
0.5 whose threshold classifications depend on floating-point roundoff. These
differences are documented in the reproduction review; they do not reverse
the primary failed practical-success decision.
