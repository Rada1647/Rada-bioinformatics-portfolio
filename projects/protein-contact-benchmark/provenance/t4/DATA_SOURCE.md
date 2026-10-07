# T4 lysozyme source and extraction

Source: https://zenodo.org/records/3989057 (version 1; DOI 10.5281/zenodo.3989057).
Original study: Kümmerer et al. (2021), *Fitting Side-Chain NMR Relaxation Data Using Molecular Simulations*, Journal of Chemical Theory and Computation 17(8), 5262–5275. DOI 10.1021/acs.jctc.0c01338. The source deposit references the earlier preprint, DOI 10.1101/2020.08.18.256024.

The deposit describes five independent 1 microsecond simulations, AMBER ff99SB*-ILDN with modified methyl rotation barriers, TIP4P/2005 water, and protein coordinates written every 1 ps. The five source trajectory names are md1us4.xtc through md1us8.xtc, with matching TPR topology files. Source-provided MD5 values and exact source sizes are preserved in zenodo_record.json.

All five TPR protein topologies have been checked to have identical CA indices, residue identities, residue IDs, and bonds. The full TPR contains 51,932 atoms; the first 2,612 atoms are the protein, with exactly 162 CA atoms. The protein-only XTC uses this contiguous protein subset. The observed 162-residue sequence is recorded in sequence_audit.json and compared with the current UniProt P00720 reference. Its Thr54 and Ala97 are consistent with the common cysteine-free T4 lysozyme construct, but the exact recorded sequence must be used rather than assuming wild type. Relative to current 164-residue P00720, other differences are R12G, I137R, and absent terminal NL; strain/reference differences should not automatically be called engineered mutations.

## Reproduction

Run acquire_extract.py with Python containing MDAnalysis 2.9.0 and NumPy 2.3.5. It processes files in numerical run order. Downloads larger than 100 MB use bounded HTTP byte ranges (initially 16 concurrent requests; then increased once to 32, capped there) from the same published Zenodo file, verifying each HTTP Content-Range header. After reconstructing each complete original file, its full MD5 must match the published checksum; the script also computes SHA256. Completion manifests retain both checksums and source size. HTTP 429 or 503 responses cause the concurrent downloader to stop, allowing a lower-concurrency restart; other range failures are logged and retried up to four attempts. Partial range downloads are resumable.

After verification, the TPR protein subset is merged into a protein-only MDAnalysis universe and paired with its corresponding XTC. Each selected frame is reconstructed using the retained covalent backbone bonds and the instantaneous periodic box, with MDAnalysis.lib.mdamath.make_whole applied to the connected N–CA–C backbone. Reconstructing unused sidechains and hydrogens is avoided. If all raw protein covalent bonds are shorter than 3 angstroms and the smallest singular value of the triclinic box matrix exceeds 6 angstroms, the molecule is already whole and reconstruction is skipped. The CA coordinates from this optimization were compared with the original full-protein reconstruction for all 5,001 frames of md1us4: all 2,430,486 coordinate scalars and all other NPZ fields were bitwise identical (unwrap_equivalence.json). The original full-protein coordinate output is retained as md1us4_full_make_whole_reference.npz. CA coordinates are sampled from source frame 0 at an exact 200 ps stride; source and sampled timestamps, frame indices, box dimensions, CA indices, residue IDs, and residue names are saved. Distances are in angstroms. No contact features, denoising estimates, or benchmark outcomes are calculated during acquisition/extraction.

The script checks finite coordinates, exact 200 ps timestamp differences, one connected protein fragment, and consecutive CA distances below 5 angstroms after reconstruction. It records the maximum CA-neighbor separation and magnitude of any periodic-coordinate correction. A completed run is indicated by its *_extraction.json manifest and matching *_ca_200ps.npz output, not by the existence of an incomplete raw download.

Full source XTC files are approximately 10 GB each, so newly downloaded raw files md1us4 through md1us7 are removed after verified successful extraction to stay within available disk. The final md1us8 raw file is retained for independent source-level review. Checksums, TPRs, source metadata, extraction code, and sampled coordinates are retained. Each raw file can be downloaded again from its recorded source. No pre-existing user data are deleted.

## License provenance

The source record marks the dataset open access but its fetched metadata supplies no explicit data-license field, and the web page's License heading does not name a license. We do not assign a license to the source data or infer a license from the separate ABSURDer software repository. Public distribution of this project should provide source download/extraction instructions and verify data redistribution terms before bundling original or derived coordinates.

## Implementation audit trail

The initial full-protein md1us4 coordinate extraction completed, but its metadata write raised TypeError because a file-size property was mistakenly called as a function. The coordinates and checksum-verified raw input were preserved. The property call was corrected before any contact benchmark; re-extraction using the optimized backbone route produced bitwise identical coordinates and successfully wrote the manifest. The retained acquisition log records this resolved failure rather than hiding it.

## Runnable commands after unpacking

From this `data/t4` directory, with the documented Python dependencies installed:

```bash
python acquire_extract.py --workers 32
```

The script locates `zenodo_record.json` relative to itself, so no original workspace paths are needed. Completed runs with extraction manifests are skipped. To independently download and re-extract all original data into a clean sibling directory, preserving the packaged outputs:

```bash
mkdir -p ../t4_reproduction
cp acquire_extract.py zenodo_record.json ../t4_reproduction/
python ../t4_reproduction/acquire_extract.py --workers 32
```

`--workers 16` or a smaller value is available for a more limited connection. The download requires approximately 50.7 GB of network transfer in total and roughly one 10.2 GB raw trajectory of working disk at a time, plus topology, index, and output files. The last raw file is retained for independent review. To extract just one specified run in a clean output directory, append for example `--runs md1us8`.

A later disk-full interruption occurred during md1us6 download. The downloader was stopped; its atomic range checkpoint and partial file were retained and checked structurally. Explicit MDAnalysis trajectory closure was added before raw-file removal, and obsolete own-task temporary files were removed separately by the coordinating process. Acquisition resumed from recorded completed spans; no partial input was accepted for extraction. `storage_interruption.json` records the event. Full-file published MD5 verification remains mandatory after resumption.
