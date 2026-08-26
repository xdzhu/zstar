"""Collect ZStar workspaces into an auditable Born-charge database."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any, Iterable

import numpy as np


SCHEMA_VERSION = "1.0"
MANIFEST_COLUMNS = (
    "material_id",
    "formula",
    "dimensionality",
    "workspace",
    "backend",
    "structure_source",
    "notes",
)


@dataclass(frozen=True)
class ManifestEntry:
    material_id: str
    formula: str
    dimensionality: int
    workspace: Path
    backend: str = "abacus-pyatb"
    structure_source: str = ""
    notes: str = ""


def _numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        clean = raw.split("#", 1)[0].strip()
        if not clean:
            continue
        fields = clean.split()
        try:
            values = [float(value.replace("D", "E").replace("d", "e")) for value in fields]
        except ValueError:
            continue
        if len(values) >= 9:
            rows.append(values[-9:])
    return rows


def read_born(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a Phonopy-style BORN file as epsilon infinity and atom tensors."""

    source = Path(path)
    rows = _numeric_rows(source)
    if len(rows) < 2:
        raise ValueError(f"BORN must contain a dielectric row and at least one tensor: {source}")
    epsilon = np.asarray(rows[0], dtype=float).reshape(3, 3)
    tensors = np.asarray(rows[1:], dtype=float).reshape(-1, 3, 3)
    return epsilon, tensors


def read_zborn(path: str | Path) -> np.ndarray:
    rows: list[list[float]] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        values: list[float] = []
        for field in raw.split():
            try:
                values.append(float(field.replace("D", "E").replace("d", "e")))
            except ValueError:
                continue
        # Z-BORN-symm.out prefixes each tensor with an atom index and symbol.
        if len(values) >= 10:
            rows.append(values[-9:])
    if not rows:
        raise ValueError(f"No 3x3 Born tensors found in {path}")
    return np.asarray(rows, dtype=float).reshape(-1, 3, 3)


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _find_first(root: Path, candidates: Iterable[str]) -> Path | None:
    for relative in candidates:
        path = root / relative
        if path.is_file():
            return path
    return None


def _read_static_response(root: Path) -> tuple[np.ndarray | None, str | None]:
    candidates = [
        "dielectric_response/ir_response_real.dat",
        "ir/ir_response_real.dat",
        "ir_response_real.dat",
        "reference/ir_response_real.dat",
    ]
    path = _find_first(root, candidates)
    if path is None:
        found = sorted(root.glob("**/ir_response_real.dat"))
        path = found[0] if found else None
    if path is None:
        return None, None
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        try:
            values = [float(item) for item in raw.split()]
        except ValueError:
            continue
        if len(values) >= 10:
            return np.asarray(values[1:10], dtype=float).reshape(3, 3), str(path)
    return None, str(path)


def _read_gap(root: Path) -> tuple[float | None, bool | None, str | None]:
    candidates = [
        root / "0.no-move" / "zstar_insulation.json",
        root / "zstar_insulation.json",
    ]
    candidates.extend(sorted(root.glob("**/zstar_insulation.json")))
    for path in candidates:
        data = _read_json(path)
        if not data:
            continue
        gap = data.get("gap_eV")
        insulating = data.get("insulating")
        return (
            float(gap) if gap is not None else None,
            bool(insulating) if insulating is not None else None,
            str(path),
        )
    return None, None, None


def _atomic_labels_from_stru(path: Path) -> list[str]:
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip().upper() == "ATOMIC_POSITIONS")
    except StopIteration:
        return []
    labels: list[str] = []
    index = start + 2
    while index < len(lines):
        symbol = lines[index].strip().split()
        if not symbol:
            index += 1
            continue
        element = symbol[0]
        if index + 2 >= len(lines):
            break
        try:
            natoms = int(lines[index + 2].split()[0])
        except (ValueError, IndexError):
            index += 1
            continue
        labels.extend([element] * natoms)
        index += 3 + natoms
    return labels


def load_manifest(path: str | Path) -> list[ManifestEntry]:
    manifest = Path(path).resolve()
    entries: list[ManifestEntry] = []
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = {"material_id", "dimensionality", "workspace"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest missing columns: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            material_id = (row.get("material_id") or "").strip()
            if not material_id:
                raise ValueError(f"Empty material_id at manifest row {row_number}")
            workspace_value = (row.get("workspace") or "").strip()
            workspace = Path(workspace_value)
            if not workspace.is_absolute():
                workspace = (manifest.parent / workspace).resolve()
            entries.append(
                ManifestEntry(
                    material_id=material_id,
                    formula=(row.get("formula") or material_id).strip(),
                    dimensionality=int(row["dimensionality"]),
                    workspace=workspace,
                    backend=(row.get("backend") or "abacus-pyatb").strip(),
                    structure_source=(row.get("structure_source") or "").strip(),
                    notes=(row.get("notes") or "").strip(),
                )
            )
    duplicate = sorted({item.material_id for item in entries if sum(x.material_id == item.material_id for x in entries) > 1})
    if duplicate:
        raise ValueError(f"Duplicate material_id values: {', '.join(duplicate)}")
    return entries


def write_manifest_template(path: str | Path) -> Path:
    target = Path(path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(MANIFEST_COLUMNS)
        writer.writerow(["bto-001", "BaTiO3", 3, "cases/3d_bulk/BaTiO3/work", "abacus-pyatb", "doi-or-database-id", "candidate"])
    return target


def collect_entry(entry: ManifestEntry) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    root = entry.workspace
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "material_id": entry.material_id,
        "formula": entry.formula,
        "dimensionality": entry.dimensionality,
        "backend": entry.backend,
        "workspace": str(root),
        "structure_source": entry.structure_source,
        "notes": entry.notes,
        "status": "incomplete",
        "quality_flags": [],
    }
    tensors: np.ndarray | None = None
    epsilon_inf: np.ndarray | None = None
    born_path = _find_first(root, ["BORN", "BORN-for-phonopy.out"])
    zborn_path = _find_first(root, ["Z-BORN-symm.out", "Z-BORN-all.out"])
    tensor_path: Path | None = zborn_path or born_path
    if born_path:
        epsilon_inf, born_tensors = read_born(born_path)
        tensors = born_tensors
    if zborn_path:
        tensors = read_zborn(zborn_path)

    stru = _find_first(root, ["STRU", "0.no-move/STRU"])
    labels = _atomic_labels_from_stru(stru or Path())
    natoms_structure = len(labels) or None
    record["natoms_structure"] = natoms_structure
    if entry.dimensionality == 3:
        record["response_kind"] = "bulk_3d"
    elif entry.dimensionality == 2:
        record["response_kind"] = "sheet_2d"
    else:
        record["response_kind"] = "molecular_spectroscopy"

    gap, insulating, gap_source = _read_gap(root)
    static_total, static_source = _read_static_response(root)
    if epsilon_inf is not None:
        record["epsilon_infinity"] = epsilon_inf.tolist()
        record["k_electronic_mean"] = float(np.trace(epsilon_inf) / 3.0)
    if static_total is not None and entry.dimensionality == 3:
        record["epsilon_static_total"] = static_total.tolist()
        record["k_static_mean"] = float(np.trace(static_total) / 3.0)
        record["high_k_rank_basis"] = "total_static_3d"
    elif epsilon_inf is not None and entry.dimensionality == 3:
        record["high_k_rank_basis"] = "electronic_only_not_ranked"
        record["quality_flags"].append("missing_total_static_dielectric")
    else:
        record["high_k_rank_basis"] = "not_applicable"

    record.update({"gap_eV": gap, "insulating": insulating})
    if insulating is False:
        record["quality_flags"].append("metallic_reference")
    if gap is None and entry.dimensionality in {2, 3}:
        record["quality_flags"].append("missing_gap_gate")

    atom_rows: list[dict[str, Any]] = []
    if tensors is not None:
        record["natoms_bec"] = int(len(tensors))
        full_cell = natoms_structure is None or len(tensors) == natoms_structure
        record["tensor_scope"] = "full_cell" if full_cell else "symmetry_representatives"
        record["max_bec_component_abs_e"] = float(np.max(np.abs(tensors)))
        record["max_bec_singular_value_e"] = float(
            max(np.linalg.svd(tensor, compute_uv=False)[0] for tensor in tensors)
        )
        if full_cell:
            acoustic = np.sum(tensors, axis=0)
            record["acoustic_sum_tensor"] = acoustic.tolist()
            record["acoustic_sum_max_abs_e"] = float(np.max(np.abs(acoustic)))
            if record["acoustic_sum_max_abs_e"] > 0.1:
                record["quality_flags"].append("large_acoustic_sum_residual")
        else:
            record["quality_flags"].append("representative_tensors_only")
        if len(labels) != len(tensors):
            labels = [f"atom-{index + 1}" for index in range(len(tensors))]
        for atom_index, (label, tensor) in enumerate(zip(labels, tensors), start=1):
            atom_rows.append(
                {
                    "schema_version": SCHEMA_VERSION,
                    "material_id": entry.material_id,
                    "atom_index": atom_index,
                    "element": label,
                    "tensor_e": tensor.tolist(),
                    "trace_over_3_e": float(np.trace(tensor) / 3.0),
                    "frobenius_norm_e": float(np.linalg.norm(tensor)),
                }
            )
        record["status"] = "complete" if insulating is not False else "rejected_metal"
        if not full_cell and record["status"] == "complete":
            record["status"] = "incomplete"
    else:
        has_spectra = any(root.glob("**/ir_modes.csv")) or any(root.glob("**/raman_modes.csv"))
        if entry.dimensionality == 0 and has_spectra:
            record["status"] = "complete_auxiliary"
            record["tensor_scope"] = "not_applicable"
        else:
            record["quality_flags"].append("missing_bec_tensor")

    record["provenance"] = {
        "born": str(tensor_path) if tensor_path else None,
        "gap": gap_source,
        "static_response": static_source,
    }
    return record, atom_rows


def collect_database(manifest: str | Path, output: str | Path) -> dict[str, Any]:
    entries = load_manifest(manifest)
    outdir = Path(output).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    tensors: list[dict[str, Any]] = []
    for entry in entries:
        record, rows = collect_entry(entry)
        records.append(record)
        tensors.extend(rows)

    columns = [
        "material_id", "formula", "dimensionality", "backend", "status",
        "gap_eV", "insulating", "natoms_structure", "natoms_bec", "tensor_scope",
        "response_kind", "k_electronic_mean",
        "k_static_mean", "high_k_rank_basis", "max_bec_component_abs_e",
        "max_bec_singular_value_e", "acoustic_sum_max_abs_e", "quality_flags",
        "workspace", "structure_source", "notes",
    ]
    with (outdir / "materials.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            flat = dict(record)
            flat["quality_flags"] = ";".join(record["quality_flags"])
            writer.writerow(flat)
    with (outdir / "materials.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
    with (outdir / "born_tensors.jsonl").open("w", encoding="utf-8") as handle:
        for row in tensors:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")

    ranked = sorted(
        (
            record for record in records
            if record.get("high_k_rank_basis") == "total_static_3d"
            and record.get("status") == "complete"
        ),
        key=lambda item: float(item["k_static_mean"]),
        reverse=True,
    )
    with (outdir / "high_k_rank.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "material_id", "formula", "k_static_mean", "gap_eV", "max_bec_singular_value_e"])
        for rank, record in enumerate(ranked, start=1):
            writer.writerow([rank, record["material_id"], record["formula"], record["k_static_mean"], record.get("gap_eV"), record.get("max_bec_singular_value_e")])

    summary = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "manifest": str(Path(manifest).resolve()),
        "materials": len(records),
        "complete": sum(record["status"] == "complete" for record in records),
        "complete_auxiliary": sum(record["status"] == "complete_auxiliary" for record in records),
        "rejected_metal": sum(record["status"] == "rejected_metal" for record in records),
        "incomplete": sum(record["status"] == "incomplete" for record in records),
        "ranked_high_k_3d": len(ranked),
        "atom_tensors": len(tensors),
        "files": ["materials.csv", "materials.jsonl", "born_tensors.jsonl", "high_k_rank.csv"],
    }
    (outdir / "database_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary
