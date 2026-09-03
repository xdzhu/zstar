"""Resolve and stage ABACUS pseudopotential and orbital assets.

The resolver deliberately never edits the user's source STRU and never uses
symlinks.  Generated displacement folders receive ordinary file copies via
the existing input-set staging code.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable


SECTION_HEADERS = {
    "ATOMIC_SPECIES",
    "NUMERICAL_ORBITAL",
    "NUMERICAL_DESCRIPTOR",
    "ABFS_ORBITAL",
    "LATTICE_CONSTANT",
    "LATTICE_VECTORS",
    "LATTICE_PARAMETERS",
    "ATOMIC_POSITIONS",
}


class AbacusAssetError(ValueError):
    """Raised when an ABACUS asset cannot be resolved unambiguously."""


@dataclass(frozen=True)
class AssetRecord:
    kind: str
    element: str
    source: Path
    name_in_stru: str

    def as_dict(self) -> dict[str, str]:
        return {
            "kind": self.kind,
            "element": self.element,
            "source": str(self.source),
            "name_in_stru": self.name_in_stru,
            "sha256": sha256_file(self.source),
        }


@dataclass(frozen=True)
class PreparedSTRU:
    path: Path
    assets: tuple[Path, ...]
    records: tuple[AssetRecord, ...]
    changed: bool


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _section_ranges(lines: list[str]) -> dict[str, tuple[int, int]]:
    starts: list[tuple[int, str]] = []
    for index, raw in enumerate(lines):
        header = raw.split("#", 1)[0].strip().upper()
        if header in SECTION_HEADERS:
            starts.append((index, header))
    ranges: dict[str, tuple[int, int]] = {}
    for position, (start, name) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        ranges[name] = (start, end)
    return ranges


def _content_indices(lines: list[str], start: int, end: int) -> list[int]:
    indices = []
    for index in range(start + 1, end):
        stripped = lines[index].split("#", 1)[0].strip()
        if stripped and not stripped.startswith("//"):
            indices.append(index)
    return indices


def _species(lines: list[str], ranges: dict[str, tuple[int, int]]) -> list[tuple[str, str, str, int]]:
    if "ATOMIC_SPECIES" not in ranges:
        raise AbacusAssetError(
            "STRU does not contain an ATOMIC_SPECIES section; "
            "ABACUS pseudopotentials cannot be resolved."
        )
    start, end = ranges["ATOMIC_SPECIES"]
    result = []
    for index in _content_indices(lines, start, end):
        tokens = lines[index].split("#", 1)[0].split()
        if len(tokens) < 3:
            raise AbacusAssetError(
                f"Invalid ATOMIC_SPECIES line {index + 1}: {lines[index].rstrip()!r}. "
                "Expected: ELEMENT MASS PSEUDOPOTENTIAL."
            )
        result.append((tokens[0], tokens[1], tokens[2], index))
    if not result:
        raise AbacusAssetError("ATOMIC_SPECIES is present but contains no species.")
    return result


def _orbital_entries(
    lines: list[str], ranges: dict[str, tuple[int, int]], count: int
) -> list[tuple[str, int]]:
    if "NUMERICAL_ORBITAL" not in ranges:
        return []
    start, end = ranges["NUMERICAL_ORBITAL"]
    indices = _content_indices(lines, start, end)
    if len(indices) != count:
        raise AbacusAssetError(
            "NUMERICAL_ORBITAL does not contain one orbital entry per atomic species "
            f"(found {len(indices)}, expected {count}). Check the STRU ordering."
        )
    return [(lines[index].split("#", 1)[0].split()[0], index) for index in indices]


def _normalise_dir(value: str | Path | None) -> Path | None:
    if value is None or not str(value).strip():
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise AbacusAssetError(f"ABACUS asset directory does not exist: {path}")
    return path


def _matching_files(directory: Path, element: str, extension: str) -> list[Path]:
    prefix = element.lower()
    result = []
    for path in directory.rglob("*"):
        if not path.is_file() or path.suffix.lower() != extension:
            continue
        stem = path.stem.lower()
        if stem == prefix or stem.startswith(prefix + "_") or stem.startswith(prefix + "."):
            result.append(path.resolve())
    return sorted(set(result), key=lambda item: str(item).lower())


def _format_candidates(candidates: Iterable[Path]) -> str:
    return "\n".join(f"    - {path}" for path in candidates)


def _resolve_one(
    *,
    kind: str,
    element: str,
    reference: str,
    stru_dir: Path,
    search_dir: Path | None,
    extension: str,
) -> Path:
    reference_path = Path(reference).expanduser()
    source_candidate = (
        reference_path.resolve()
        if reference_path.is_absolute()
        else (stru_dir / reference_path).resolve()
    )

    # An explicitly supplied directory takes precedence for relative STRU
    # references.  Exact basename matching makes common multi-version libraries
    # deterministic without guessing from a prefix.
    if search_dir is not None:
        exact_candidates = [
            (search_dir / reference_path).resolve(),
            (search_dir / reference_path.name).resolve(),
        ]
        for candidate in exact_candidates:
            if candidate.is_file() and candidate.suffix.lower() == extension:
                return candidate
        candidates = _matching_files(search_dir, element, extension)
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise AbacusAssetError(
                f"Multiple {kind} files match element {element!r} in {search_dir}:\n"
                f"{_format_candidates(candidates)}\n"
                "Specify the exact filename in STRU, or point "
                f"--{'pp' if kind == 'pseudopotential' else 'orb'} to a narrower directory."
            )

    if source_candidate.is_file():
        return source_candidate

    option = "--pp" if kind == "pseudopotential" else "--orb"
    config_key = "abacus.pseudo_dir" if kind == "pseudopotential" else "abacus.orbital_dir"
    searched = f" under {search_dir}" if search_dir is not None else ""
    raise AbacusAssetError(
        f"Cannot find {kind} for element {element!r}: {reference!r}{searched}.\n"
        f"Use {option} DIRECTORY, set {config_key} in the global ZStar config, "
        "or correct the filename in STRU."
    )


def _replace_species_line(line: str, filename: str) -> str:
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    comment = ""
    if "#" in body:
        body, comment = body.split("#", 1)
        comment = "  #" + comment
    tokens = body.split()
    tokens[2] = filename
    return " ".join(tokens) + comment + newline


def _replace_orbital_line(line: str, filename: str) -> str:
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    comment = ""
    if "#" in body:
        body, comment = body.split("#", 1)
        comment = "  #" + comment
    leading = re.match(r"^\s*", body).group(0)
    return leading + filename + comment + newline


def prepare_stru_assets(
    stru: str | Path,
    *,
    pp_dir: str | Path | None = None,
    orb_dir: str | Path | None = None,
    output_dir: str | Path = ".zstar",
) -> PreparedSTRU:
    """Resolve assets and write a generated STRU copy when needed.

    Relative paths in the source STRU are interpreted relative to the source
    STRU directory.  Explicit directories override relative references; an
    exact filename is preferred, followed by a unique element-prefix match.
    """
    source = Path(stru).expanduser().resolve()
    pp_root = _normalise_dir(pp_dir)
    orb_root = _normalise_dir(orb_dir)
    if pp_root is None and orb_root is None:
        # Preserve the legacy route when the user did not request asset
        # resolution.  Existing ABACUS examples may intentionally rely on
        # their own runtime environment.
        return PreparedSTRU(source, tuple(), tuple(), False)
    if not source.is_file():
        raise AbacusAssetError(f"STRU file does not exist: {source}")
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    ranges = _section_ranges(lines)
    species = _species(lines, ranges)
    orbitals = _orbital_entries(lines, ranges, len(species))

    records: list[AssetRecord] = []
    replacements: dict[int, str] = {}
    assets: list[Path] = []
    for element, _mass, reference, index in species:
        path = _resolve_one(
            kind="pseudopotential",
            element=element,
            reference=reference,
            stru_dir=source.parent,
            search_dir=pp_root,
            extension=".upf",
        )
        records.append(AssetRecord("pseudopotential", element, path, path.name))
        replacements[index] = _replace_species_line(lines[index], path.name)
        if path not in assets:
            assets.append(path)

    for (element, _mass, _reference, _index), (reference, index) in zip(species, orbitals):
        path = _resolve_one(
            kind="orbital",
            element=element,
            reference=reference,
            stru_dir=source.parent,
            search_dir=orb_root,
            extension=".orb",
        )
        records.append(AssetRecord("orbital", element, path, path.name))
        replacements[index] = _replace_orbital_line(lines[index], path.name)
        if path not in assets:
            assets.append(path)

    changed_lines = list(lines)
    for index, replacement in replacements.items():
        changed_lines[index] = replacement
    changed = changed_lines != lines
    destination_dir = Path(output_dir).expanduser()
    if not destination_dir.is_absolute():
        destination_dir = source.parent / destination_dir
    destination_dir = destination_dir.resolve()
    destination_dir.mkdir(parents=True, exist_ok=True)
    if changed:
        destination = destination_dir / "STRU.resolved"
        destination.write_text("".join(changed_lines), encoding="utf-8", newline="\n")
    else:
        destination = source

    manifest = {
        "source_stru": str(source),
        "prepared_stru": str(destination),
        "changed": changed,
        "assets": [record.as_dict() for record in records],
    }
    (destination_dir / "assets.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return PreparedSTRU(destination, tuple(assets), tuple(records), changed)
