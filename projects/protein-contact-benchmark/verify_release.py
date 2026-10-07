"""Check packaged file hashes without installing scientific dependencies."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def main():
    manifest = json.loads((ROOT / "RELEASE_MANIFEST.json").read_text())
    errors = []
    for record in manifest["files"]:
        relative = Path(record["path"])
        if relative.is_absolute() or ".." in relative.parts:
            errors.append({"path": str(relative), "issue": "unsafe manifest path"})
            continue
        path = ROOT / relative
        if not path.is_file():
            errors.append({"path": str(relative), "issue": "missing"})
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if path.stat().st_size != record["bytes"] or digest != record["sha256"]:
            errors.append({"path": str(relative), "issue": "changed"})
    print(json.dumps({"files_checked": len(manifest["files"]), "errors": errors,
                      "status": "pass" if not errors else "fail"}, indent=2))
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
