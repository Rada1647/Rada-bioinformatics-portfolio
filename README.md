# Rada's Bioinformatics Portfolio

Computational biology projects exploring protein structure, function, and the
methods used to analyze biological data. My interests include proteins and
antibodies, and this portfolio documents my progress through practical research
questions.

Each completed project records the question, data sources, methods, results,
limitations, and steps needed to check the work. Negative findings are included
alongside any future positive results.

## Projects

| Project | Question | Status and result |
|---|---|---|
| [Protein-contact smoothing benchmark](projects/protein-contact-benchmark/) | Can an expected-median mathematical construction improve recovery of artificially corrupted protein contacts? | Completed computational benchmark. The predefined practical-success criterion failed on T4 lysozyme. A separate villin stress test also favored stronger controls. |

### Protein-contact smoothing benchmark

The first project tests a construction discussed in OpenAI/math catalogue family
#332 against conventional smoothing and sharpening controls. It uses deposited
molecular-dynamics trajectories, a privately frozen protocol, and held-out data.

The cubic method's mean Brier error was **8.53% higher** than the validation-selected
control for the primary T4 experiment and **24.79% higher** in the separate villin
stress test. Ubiquitin was excluded because the available temporal provenance
did not meet the declared requirements. This result concerns the tested engineering
application; it does not disprove the source mathematical theorem.

- [Project overview and results](projects/protein-contact-benchmark/README.md)
- [Research report (PDF)](projects/protein-contact-benchmark/benchmark_report.pdf)
- [Reproduction instructions](projects/protein-contact-benchmark/REPRODUCING.md)
- [Frozen rules](projects/protein-contact-benchmark/protocol/FROZEN_RULES.txt)
- [Data sources and attribution](projects/protein-contact-benchmark/DATA_SOURCES.md)
- [References](projects/protein-contact-benchmark/REFERENCES.txt)

## Run the lightweight checks

From a clone of this repository, using Python 3.12:

```bash
cd projects/protein-contact-benchmark
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-benchmark.txt
.venv/bin/python scripts/check_release.py
```

These offline checks verify file hashes, small mathematical examples, synthetic
benchmark behavior, and regeneration of the archived summary. They do not
download source datasets or rerun the full protein experiment. Full reconstruction
requires separate source downloads; the project guide explains the procedure
and resource requirements.

## Scope and transparency

This is a personal learning and research portfolio. The benchmark was developed
with AI assistance and checked through five independent AI technical review roles.
Those checks are documented in the [review records](projects/protein-contact-benchmark/reviews/);
they are not human peer review or experimental biological validation.

Only the completed project above is presented as a result. Future studies and
learning exercises can be added in their own folders with clear status labels.
Source trajectories, topologies, and coordinate/contact arrays are not included.

Original project code and writing use the [MIT license](LICENSE). Third-party
datasets, metadata, papers, and dependencies retain their own terms; see the
project's [notices](projects/protein-contact-benchmark/THIRD_PARTY_NOTICES.md).
