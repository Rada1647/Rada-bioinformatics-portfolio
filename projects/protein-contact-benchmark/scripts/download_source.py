#!/usr/bin/env python3
"""Acquire one pinned archive; never treat a partial or unchecked file as data."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "villin": {
        "name": "protein_folding_datasets.zip",
        "url": "https://ndownloader.figshare.com/files/65180772",
        "bytes": 2559410654,
        "md5": "484eafab71b27e3d66f48d4820d08f47",
        "sha256": "4f512459ac7bfc42e8bcd6b247648b9ccdfe60dde7d8ace1dbeb765293618fff",
        "source": "https://doi.org/10.6084/m9.figshare.32541291.v1",
    },
    "ubiquitin": {
        "name": "ubiquitin-md-generated-ensemble.zip",
        "url": "https://zenodo.org/records/7792288/files/ubiquitin-md-generated-ensemble.zip?download=1",
        "bytes": 228792237,
        "md5": "3372f857beb371ca5ebc8d57877c5c79",
        "sha256": "d21f88be709745f858a78c11dbfb7ff6bdb014e6cb0aaac2b489071f540dd0d8",
        "source": "https://doi.org/10.5281/zenodo.7792288",
    },
}


def hashes(path: Path) -> dict[str, str]:
    hs = {name: hashlib.new(name) for name in ("md5", "sha256")}
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            for digest in hs.values():
                digest.update(block)
    return {name: digest.hexdigest() for name, digest in hs.items()}


def verify(path: Path, source: dict) -> dict[str, str]:
    if path.stat().st_size != source["bytes"]:
        raise RuntimeError(f"Unexpected archive size for {path}; do not analyze this file")
    result = hashes(path)
    if any(result[name] != source[name] for name in ("md5", "sha256")):
        raise RuntimeError(f"Archive checksum mismatch for {path}; do not analyze this file")
    return result


def acquire(dataset: str) -> Path:
    source = SOURCES[dataset]
    directory = ROOT / "data" / dataset
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / source["name"]
    if target.exists():
        result = verify(target, source)
        print(f"Verified existing {target}", flush=True)
    else:
        partial = target.with_suffix(target.suffix + ".part")
        # Failed downloads can be inspected or removed explicitly. Never append
        # a new response to an unverified partial file.
        if partial.exists():
            raise RuntimeError(f"Incomplete download exists: {partial}. Inspect/remove it before retrying.")
        request = urllib.request.Request(source["url"], headers={"User-Agent": "protein-contact-negative-benchmark/1.0"})
        print(f"Downloading {source['bytes']:,} bytes from {source['source']}", flush=True)
        with urllib.request.urlopen(request, timeout=180) as response, partial.open("xb") as output:
            for block in iter(lambda: response.read(8 * 1024**2), b""):
                output.write(block)
        result = verify(partial, source)
        partial.replace(target)
    receipt = {
        **source,
        "verified_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verified_hashes": result,
        "purpose": "Local reproduction only; no source archive is part of the repository",
    }
    (directory / "download_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=tuple(SOURCES))
    args = parser.parse_args()
    acquire(args.dataset)


if __name__ == "__main__":
    main()
