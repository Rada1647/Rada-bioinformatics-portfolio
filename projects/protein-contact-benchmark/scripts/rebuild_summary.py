#!/usr/bin/env python3
"""Rebuild archived summary/CSV in isolation and compare their scientific content.

No coordinates, source downloads, or live data/prepared directory are needed.
This checks report aggregation, not the correctness of the underlying simulation
or benchmark scores. The frozen collector executes unchanged in a temporary tree.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="Optional new directory in which to retain regenerated summary.json and metrics.csv")
    args = parser.parse_args()
    if args.out is not None and args.out.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {args.out}")
    with tempfile.TemporaryDirectory(prefix="protein332_summary_") as temporary:
        work = Path(temporary)
        paths = ["collect_results.py", "protocol/protocol_v1.1.json"]
        for dataset in ("t4", "villin"):
            paths.extend(str(p.relative_to(ROOT)) for p in sorted((ROOT / "results" / dataset).glob("*.json")))
        for relative in paths:
            target = work / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative, target)
        for dataset in ("t4", "villin"):
            for name in ("manifest.json", "coverage.json"):
                target = work / "data" / "prepared" / dataset / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / "provenance" / "prepared" / dataset / name, target)
        qualification = work / "data" / "ubiquitin" / "qualification.json"
        qualification.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "provenance" / "ubiquitin" / "qualification.json", qualification)
        completed = subprocess.run([sys.executable, str(work / "collect_results.py")], cwd=work, text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError(f"Frozen collector failed:\n{completed.stdout}\n{completed.stderr}")
        expected = json.loads((ROOT / "summary.json").read_text())
        regenerated = json.loads((work / "summary.json").read_text())
        # This collector writes exactly one fresh aggregation timestamp.
        expected.pop("generated_utc", None)
        regenerated.pop("generated_utc", None)
        summary_match = expected == regenerated
        csv_match = (ROOT / "metrics.csv").read_bytes() == (work / "metrics.csv").read_bytes()
        if args.out is not None:
            args.out.mkdir(parents=True, exist_ok=False)
            for name in ("summary.json", "metrics.csv"):
                shutil.copyfile(work / name, args.out / name)
        print(json.dumps({
            "summary_scientific_content_exact_match": summary_match,
            "metrics_csv_byte_exact_match": csv_match,
            "excluded_summary_fields": ["generated_utc"],
            "scope": "Reaggregation of archived numerical outputs, not a benchmark rerun",
        }, indent=2))
        if not summary_match or not csv_match:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
