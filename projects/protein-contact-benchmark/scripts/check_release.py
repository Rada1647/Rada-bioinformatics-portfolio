#!/usr/bin/env python3
"""Run the offline release checks; no protein data or network access is used."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    commands = [
        ["verify_release.py"],
        ["check_methods.py"],
        ["benchmark_extension.py", "--self-test"],
        ["scripts/rebuild_summary.py"],
    ]
    for command in commands:
        print("Running " + " ".join(command), flush=True)
        subprocess.run([sys.executable, *command], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
