#!/usr/bin/env python3
"""Reinspect the pinned ubiquitin archive's structure and DCD header metadata.

This checks archive evidence, not the historical literature/repository search.
It never invents temporal runs or computes denoising/contact benchmark scores.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import struct
import warnings
import zipfile

import MDAnalysis as mda
import numpy as np

from download_source import ROOT, acquire


def main():
    source = acquire("ubiquitin")
    out = ROOT / "data" / "ubiquitin"
    extracted = out / "archive"
    if extracted.exists():
        raise RuntimeError(f"Refusing to reuse an existing extracted directory: {extracted}. Inspect/move it before repeating.")
    with zipfile.ZipFile(source) as archive:
        members = archive.infolist()
        for member in members:
            name = PurePosixPath(member.filename)
            if name.is_absolute() or ".." in name.parts or "\\" in member.filename:
                raise RuntimeError("Unsafe ZIP member path")
            if ((member.external_attr >> 16) & 0o170000) == 0o120000:
                raise RuntimeError("ZIP symbolic links are not accepted")
        if archive.testzip() is not None:
            raise RuntimeError("ZIP CRC failure")
        archive.extractall(extracted)
    pdbs = list(extracted.rglob("topology.pdb"))
    if len(pdbs) != 1:
        raise RuntimeError("Expected one topology.pdb")
    pdb = pdbs[0]
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        universe = mda.Universe(str(pdb))
        ca = universe.select_atoms("protein and name CA")
        topology = {
            "file": str(pdb.relative_to(out)),
            "sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
            "n_atoms": universe.atoms.n_atoms,
            "n_residues": len(universe.residues),
            "n_CA": len(ca),
        }
        trajectories = []
        for dcd in sorted(extracted.rglob("*.dcd")):
            with dcd.open("rb") as stream:
                length = struct.unpack("<i", stream.read(4))[0]
                if length != 84:
                    raise RuntimeError("Unexpected DCD header record")
                header = stream.read(length)
            current = mda.Universe(str(pdb), str(dcd))
            trajectory = current.trajectory
            trajectory[0]
            first = float(trajectory.time)
            trajectory[-1]
            last = float(trajectory.time)
            trajectories.append({
                "file": str(dcd.relative_to(out)),
                "sha256": hashlib.sha256(dcd.read_bytes()).hexdigest(),
                "n_atoms": trajectory.n_atoms,
                "n_frames": len(trajectory),
                "reader_dt_ps": float(trajectory.dt),
                "reader_first_time_ps": first,
                "reader_last_time_ps": last,
                "header_nset": struct.unpack("<i", header[4:8])[0],
                "header_istart": struct.unpack("<i", header[8:12])[0],
                "header_nsavc": struct.unpack("<i", header[12:16])[0],
                "header_delta_raw_float": struct.unpack("<f", header[40:44])[0],
            })
            trajectory.close()
    result = {
        "scope": "Recomputed archive-only metadata; no literature or author-repository search is rerun",
        "source": "https://doi.org/10.5281/zenodo.7792288",
        "software": {"MDAnalysis": mda.__version__, "numpy": np.__version__},
        "archive_members": [member.filename for member in members],
        "topology": topology,
        "trajectories": trajectories,
        "warnings": sorted(set(str(item.message) for item in seen)),
        "interpretation": "File-wide DCD header timing cannot establish original run boundaries or chronological frame provenance. Consult the recorded qualification and source review before inferring temporal eligibility.",
        "denoising_metrics_computed": False,
    }
    destination = out / "archive_metadata_reproduced.json"
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(destination), "topology": topology, "trajectories": len(trajectories)}, indent=2))


if __name__ == "__main__":
    main()
