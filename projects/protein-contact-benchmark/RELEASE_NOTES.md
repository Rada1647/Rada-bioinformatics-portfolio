# Release 0.1.0 — negative-result benchmark

This benchmark is included in Rada’s Bioinformatics Portfolio. The following
record describes the original standalone preparation. Historical statements
about unpublished status refer to that earlier preparation task. The included report is the
completed research report, not a newly peer-reviewed publication. The primary
T4 practical-success criterion failed. Villin is a separate stress experiment;
ubiquitin failed temporal qualification and was not denoised.

## Preserved scientific record

The production benchmark (`benchmark_extension.py`, `methods.py`), mathematical
checks, contact preparation, scenario runner, reproduction comparator, collector,
report builder, frozen protocol, and all 34 validation/test/clean output files
are copied without changing their bytes. `provenance/RELEASE_ORIGINS.json` records
original and public paths with hashes. The archived summary, CSV, figure, report,
references and response to reviewers are also preserved.

Protocol v1.1 SHA-256:
`1e25ff2c9890809f5432675c259f8ec705d2ce9bcca23969111f771fd520c84c`.

Private prospective freezing is not public preregistration. A hash demonstrates
file identity, not an independently certified historical timestamp. Earlier
adenylate-kinase work supplied development information, including the fixed
matched-budget alpha=3 control; it is not additional prospective evidence in
this release.

## Packaging changes

- Added GitHub-facing README, reproduction instructions, attribution, citation
  metadata, MIT license for original work, contribution guidance and offline CI.
- Added verified source-download and qualification helpers plus an isolated
  summary-rebuilding helper. They do not alter the denoising methods or outcomes.
- The public T4 acquisition script checks that an extracted file actually exists
  and matches its recorded hash before accepting a completion record. Original
  extraction calculations are unchanged. The Villin extractor is byte-identical
  to the original. Source extractor changes and hashes are recorded separately.
- Source completion records and prepared-input metadata live in `provenance/`.
  They do not populate live `data/prepared/` or cause a fresh download to skip
  missing arrays.
- Removed structural source files, derived arrays, third-party article PDFs,
  downloaded repositories, duplicated reproduction outputs, temporary files,
  invocation logs and private storage receipts from this release's allowlist.
- In `reviews/03_statistics_final.json`, eight private storage-ID/metadata fields
  were replaced with explicit redaction markers. Scientific values were not
  changed. The original and redacted hashes and JSON paths are recorded in
  `provenance/RELEASE_ORIGINS.json`. Other included audit records are unchanged.

## Reproduction scope

The earlier private package bundled derived arrays. Historical references to
bundled arrays in the report or review records do not describe this public
package. Follow `REPRODUCING.md` here. Source data are acquired separately,
including about 50.7 GB for T4 and 2.56 GB for villin. Ubiquitin download is
optional and supports qualification only.

The original study's same-host, fresh-environment reruns compared 1,175,512
scientific numeric fields with zero differences. Its five AI reviews and their
limits are included. Release checks verify file identity, small mathematical
examples, synthetic integration and reaggregation of archived results. They
do not claim a new full raw-data reconstruction or independent-machine rerun.

The GitHub workflow is supplied but has not run on GitHub before publication.
Any local release checks are documented in `release_validation.json`.

## Maintaining the archive

`RELEASE_MANIFEST.json` is a snapshot of this prepared release. If packaging files
are deliberately changed, review those changes and regenerate the manifest;
never regenerate it merely to conceal an unexplained scientific-file mismatch.
Original scientific identities are additionally pinned in `RELEASE_ORIGINS.json`.
New experiments should preserve these files and use a separate protocol/version.

## Portfolio integration

The project was placed under `projects/protein-contact-benchmark/` in
`Rada1647/Rada-bioinformatics-portfolio`. The overview and citation metadata
now point to that repository. The active GitHub workflow lives in the portfolio
root and uses this project as its working directory. The nested workflow is
retained as part of the standalone package; GitHub does not run it from here.
No benchmark algorithm, frozen protocol or numerical result was changed.
The public file manifest was refreshed for these documentation changes.
