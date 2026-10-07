#!/usr/bin/env python3
"""Metadata-only qualification and C-alpha extraction; computes no contacts.

Every XTC is preserved as its own NPZ. Do not combine trajectory time axes.
Dependencies: MDAnalysis 2.9.0 and NumPy. Dataset: Figshare 32541291.v1.
"""
from pathlib import Path, PurePosixPath
from collections import Counter, defaultdict
import argparse
import hashlib
import json
import re
import tempfile
import zipfile
import warnings

import MDAnalysis as mda
import numpy as np


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, default=Path(__file__).with_name("protein_folding_datasets.zip"))
    parser.add_argument("--out", type=Path, default=Path(__file__).parent)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((out / "figshare_metadata.json").read_text())
    source = next(x for x in meta["files"] if x["name"] == args.archive.name)
    md5 = digest(args.archive, "md5")
    assert args.archive.stat().st_size == source["size"], "Archive byte count mismatch"
    assert md5 == source["computed_md5"], "Archive Figshare MD5 mismatch"
    source_sha = digest(args.archive)
    manifest = []
    topology_details = {}
    with zipfile.ZipFile(args.archive) as archive:
        infos = sorted(archive.infolist(), key=lambda x: x.filename)
        write_json(out / "archive_manifest.json", [dict(path=i.filename, size=i.file_size,
                   compressed_size=i.compress_size, crc32=f"{i.CRC:08x}") for i in infos])
        pdbs = [i.filename for i in infos if i.filename.endswith(".pdb")]
        xtcs = [i.filename for i in infos if i.filename.endswith(".xtc")]
        others = [i.filename for i in infos if not i.is_dir() and not i.filename.endswith((".pdb", ".xtc"))]
        print(json.dumps(dict(xtcs=len(xtcs), pdbs=pdbs, other_files=others)), flush=True)
        for path in pdbs:
            p = PurePosixPath(path)
            group = p.parts[1]
            raw = archive.read(path)
            dst = out / "topologies" / f"group_{group}_{p.name}"
            dst.parent.mkdir(exist_ok=True)
            dst.write_bytes(raw)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(dst))
            ca = u.select_atoms("name CA")
            info = dict(source_path=path, local_path=str(dst.relative_to(out)),
                        sha256=hashlib.sha256(raw).hexdigest(), atoms=len(u.atoms),
                        n_ca=len(ca), ca_indices=ca.indices.tolist(),
                        ca_resids=ca.resids.tolist(), ca_resnames=ca.resnames.tolist(),
                        ca_chainids=ca.chainIDs.tolist(), ca_segids=ca.segids.tolist(),
                        all_resnames=u.residues.resnames.tolist())
            assert len(set(ca.chainIDs)) == 1 and len(set(ca.segids)) == 1
            assert len(set(ca.resids)) == len(ca) == 35
            assert np.all(np.diff(ca.resids) == 1), "Missing CA residue"
            topology_details[group] = info
        write_json(out / "topology_metadata.json", topology_details)
        for path in xtcs:
            parts = PurePosixPath(path).parts
            group, folder = parts[1], parts[3]
            child = re.match(r"e(\d+)s(\d+)", folder)
            parent = re.search(r"_e(\d+)s(\d+)p(\d+)f(\d+)", folder)
            rec = dict(source_path=path, group=group, folder=folder,
                       trajectory_id=f"group{group}__{folder}",
                       child_node=f"{group}/e{child[1]}s{child[2]}" if child else None,
                       epoch=int(child[1]) if child else None,
                       sample=int(child[2]) if child else None,
                       parent_node=f"{group}/e{parent[1]}s{parent[2]}" if parent else None,
                       parent_piece=int(parent[3]) if parent else None,
                       parent_frame=int(parent[4]) if parent else None)
            manifest.append(rec)
        # The source naming convention encodes parent pointers; report missing
        # pointers explicitly, rather than inventing independent roots.
        nodes = [r["child_node"] for r in manifest]
        assert len(nodes) == len(set(nodes)), "Duplicate child keys: inspect source grouping"
        lookup = {r["child_node"]: r for r in manifest}
        for rec in manifest:
            node = rec["child_node"]
            seen = set()
            while node in lookup and lookup[node]["parent_node"] is not None:
                assert node not in seen, "Ancestry cycle"
                seen.add(node)
                node = lookup[node]["parent_node"]
            rec["ancestry_root_or_missing_parent"] = node
            rec["ancestry_complete_to_root"] = node in lookup
        write_json(out / "trajectory_manifest_metadata.json", manifest)
        summary = dict(dataset_doi=meta["doi"], dataset_license=meta["license"],
                       archive_bytes=args.archive.stat().st_size, archive_md5=md5,
                       archive_sha256=source_sha, n_trajectories=len(manifest),
                       groups=dict(Counter(r["group"] for r in manifest)),
                       epochs_by_group={g: sorted(set(r["epoch"] for r in manifest if r["group"] == g)) for g in topology_details},
                       ancestry_missing_count=sum(not r["ancestry_complete_to_root"] for r in manifest),
                       roots=dict(Counter(r["ancestry_root_or_missing_parent"] for r in manifest)),
                       software=dict(MDAnalysis=mda.__version__, numpy=np.__version__),
                       inspected_outcomes=False)
        write_json(out / "qualification_metadata.json", summary)
        print(json.dumps(summary, indent=2), flush=True)
        if args.metadata_only:
            return
        (out / "ca").mkdir(exist_ok=True)
        for k, rec in enumerate(manifest):
            topo = topology_details[rec["group"]]
            raw = archive.read(rec["source_path"])
            rec["source_xtc_sha256"] = hashlib.sha256(raw).hexdigest()
            rec["source_xtc_bytes"] = len(raw)
            with tempfile.TemporaryDirectory(dir=out, prefix="extract_") as tmp:
                xtc = Path(tmp) / "trajectory.xtc"
                xtc.write_bytes(raw)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    u = mda.Universe(str(out / topo["local_path"]), str(xtc))
                assert len(u.atoms) == topo["atoms"]
                ca = u.select_atoms("name CA")
                assert ca.indices.tolist() == topo["ca_indices"]
                times_all = np.array([ts.time for ts in u.trajectory], dtype=np.float64)
                assert np.isfinite(times_all).all()
                dt = np.diff(times_all)
                assert len(dt) and (dt > 0).all(), "Nonincreasing time axis"
                assert np.allclose(dt, 100.0, rtol=0, atol=0.1), "Expected 100 ps cadence"
                frames = np.arange(0, len(u.trajectory), 2, dtype=np.int64)
                xyz = np.empty((len(frames), len(ca), 3), dtype=np.float32)
                box = np.empty((len(frames), 6), dtype=np.float32)
                for j, frame in enumerate(frames):
                    ts = u.trajectory[int(frame)]
                    xyz[j] = ca.positions
                    box[j] = ts.dimensions if ts.dimensions is not None else np.nan
                assert np.isfinite(xyz).all(), "Nonfinite CA coordinates"
                adjacent = np.linalg.norm(np.diff(xyz.astype(np.float64), axis=1), axis=2)
                # Geometry-only sanity check for source periodic wrapping. No
                # nonlocal contact, frequency, fold, or method score is computed.
                broken_frames = np.flatnonzero(np.any(adjacent > 5.0, axis=1))
                dst = out / "ca" / (rec["trajectory_id"] + ".npz")
                np.savez_compressed(dst, xyz=xyz, time_ps=times_all[frames],
                                    source_frame=frames, dimensions_angstrom_degrees=box,
                                    ca_resids=np.asarray(topo["ca_resids"], dtype=np.int64),
                                    ca_resnames=np.asarray(topo["ca_resnames"], dtype="U8"),
                                    ca_indices=np.asarray(topo["ca_indices"], dtype=np.int64))
                rec.update(dict(npz_path=str(dst.relative_to(out)), npz_sha256=digest(dst),
                           n_source_frames=len(times_all), n_frames=len(frames), n_ca=len(ca),
                           source_time_first_ps=float(times_all[0]), source_time_last_ps=float(times_all[-1]),
                           source_dt_min_ps=float(dt.min()), source_dt_max_ps=float(dt.max()),
                           time_first_ps=float(times_all[frames[0]]), time_last_ps=float(times_all[frames[-1]]),
                           dt_ps=200.0, stride=2, coordinate_units="angstrom",
                           adjacent_ca_min_angstrom=float(adjacent.min()),
                           adjacent_ca_max_angstrom=float(adjacent.max()),
                           adjacent_ca_over_5_angstrom_frames=broken_frames.tolist(),
                           eligibility_min20_frames=len(frames) >= 20,
                           geometry_pass=not len(broken_frames)))
                u.trajectory.close()
            if k % 100 == 0:
                print(f"Extracted {k+1}/{len(manifest)} trajectories", flush=True)
        write_json(out / "trajectory_manifest.json", manifest)
        lookup = {r["child_node"]: r for r in manifest}
        ancestry_failures = []
        for rec in manifest:
            if rec["parent_node"] is not None:
                parent = lookup[rec["parent_node"]]
                if rec["parent_frame"] >= parent["n_source_frames"] or rec["epoch"] <= parent["epoch"]:
                    ancestry_failures.append(dict(child=rec["child_node"], parent=rec["parent_node"],
                                             parent_frame=rec["parent_frame"], parent_nframes=parent["n_source_frames"]))
        write_json(out / "ancestry_bounds_audit.json", dict(parent_bounds_or_epoch_failures=ancestry_failures,
                   all_pointers_present=True, n_trajectories=len(manifest),
                   n_parent_links=sum(r["parent_node"] is not None for r in manifest)))
        counts = Counter(r["n_source_frames"] for r in manifest)
        summary.update(dict(source_frame_counts=dict(counts),
                       output_frames=sum(r["n_frames"] for r in manifest),
                       output_stride=2, output_dt_ps=200, coordinate_units="angstrom",
                       trajectory_boundaries_preserved=True,
                       total_recorded_source_span_ns=sum((r["source_time_last_ps"]-r["source_time_first_ps"])/1000 for r in manifest),
                       modal_source_frames=counts.most_common(1)[0][0]))
        summary.update(dict(geometry_fail_trajectories=[r["trajectory_id"] for r in manifest if not r["geometry_pass"]],
                            fewer_than20_frames=[r["trajectory_id"] for r in manifest if not r["eligibility_min20_frames"]]))
        write_json(out / "qualification.json", summary)
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
