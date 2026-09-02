"""Small persistent manifests shared by the canonical CLI families."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Mapping


MANIFEST_SCHEMA = "zstar-workflow"
MANIFEST_VERSION = 1


def manifest_path(root: str | Path = ".", family: str = "bec") -> Path:
    return Path(root).resolve() / ".zstar" / f"{family}.json"


def write_manifest(
    family: str,
    *,
    root: str | Path = ".",
    calculator: str,
    dimensionality: int,
    options: Mapping[str, Any] | None = None,
) -> Path:
    target = manifest_path(root, family)
    target.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema": MANIFEST_SCHEMA,
        "schema_version": MANIFEST_VERSION,
        "family": str(family),
        "calculator": str(calculator).lower(),
        "dimensionality": int(dimensionality),
        "root": str(Path(root).resolve()),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "options": dict(options or {}),
    }
    target.write_text(json.dumps(data, indent=2), encoding="utf-8", newline="\n")
    return target


def read_manifest(family: str, root: str | Path = ".") -> dict[str, Any]:
    target = manifest_path(root, family)
    if not target.is_file():
        raise FileNotFoundError(
            f"No {family} workflow manifest found at {target}; run `zstar {family} pre` first"
        )
    data = json.loads(target.read_text(encoding="utf-8"))
    if data.get("schema") != MANIFEST_SCHEMA or data.get("schema_version") != MANIFEST_VERSION:
        raise ValueError(f"Unsupported ZStar workflow manifest: {target}")
    if data.get("family") != family:
        raise ValueError(f"Expected {family!r} workflow manifest: {target}")
    return data
