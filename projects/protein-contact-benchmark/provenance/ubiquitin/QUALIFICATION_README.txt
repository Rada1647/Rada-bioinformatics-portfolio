PROTHON ubiquitin qualification for protein-contact benchmark #332

Decision: NOT QUALIFIED for the confirmatory temporal benchmark.
This is a provenance limitation, not a claim that the published ensemble data are invalid.

Downloaded: https://doi.org/10.5281/zenodo.7792288 (version 1).
Archive: ubiquitin-md-generated-ensemble.zip, 228,792,237 bytes.
Published and verified MD5: 3372f857beb371ca5ebc8d57877c5c79.
SHA-256: d21f88be709745f858a78c11dbfb7ff6bdb014e6cb0aaac2b489071f540dd0d8.
Dataset license: CC-BY-4.0 according to the saved Zenodo API metadata.
Creators: Adekunle Aina, Shawn C. C. Hsueh, Steven S. Plotkin.
Paper: Aina et al., J. Chem. Inf. Model. 2023, 63(11), 3453-3461.
DOI: 10.1021/acs.jcim.3c00145.

The archive contains topology.pdb and Q75, Q80, Q85, Q90, Q95, Q99 DCDs.
The topology has 1,231 atoms and 76 residues/CA atoms with the expected
ubiquitin sequence. Each DCD has 5,000 conformations. Q99 is the native
reference ensemble; the other five are biased partially folded ensembles.
They are not six independent simulations of native ubiquitin.

The article reports five native simulation runs, each 50 ns, retaining the
final 40 ns sampled every 40 ps. However, it does not specify how those runs
map onto archived Q99 frame indices. The ZIP has no README, assembly script,
or run manifest. Original and current author repositories did not supply
such a mapping in the inspected files. A current reanalysis calls Q99 a
trajectory and uses contiguous prefixes, but that does not identify the
original independent-run boundaries. Frame-count arithmetic alone is not
sufficient evidence for five consecutive 1,000-frame segments.

The DCD headers give DELTA=1, NSAVC=1, ISTART=0, and generic DCD-plugin
titles. MDAnalysis 2.9.0 interprets the interval as 0.04888821 ps, inconsistent
with the paper's 40 ps sampling; these file-wide header constants cannot
recover original run IDs or validate physical frame times.

Consequently no runs were inferred, no per-run arrays were exported, and
no contact/denoising benchmark results were computed. To qualify the dataset,
obtain a reliable assembly script/run manifest, or the original five separate
native trajectories with confirmed frame order and timestamps.

Audit/reproduction
  python download_qualification.py
  python inspect_archive.py

Run from any directory; paths are resolved relative to the scripts. Inspection
requires NumPy and MDAnalysis (used versions: 2.3.5 and 2.9.0). Acquisition uses
only Python's standard library. The saved zenodo_record.json freezes original
metadata; the acquisition script fetches it if absent. The first full-stream
attempt was truncated, detected by checksum, and recovered with verified HTTP
byte ranges via download_ranges.py. Extraction checks member paths and ZIP CRCs.

Machine-readable decision: qualification.json.
Archive/structure metadata: download_audit.json, archive_metadata.json.
Evidence pointers and pinned author-repository commits: source_provenance.json.
The source article PDF and cloned repositories are research working copies;
they need not be included in a compact benchmark reproduction package. Include
the acquisition/inspection scripts, metadata, decision, and this README instead.
