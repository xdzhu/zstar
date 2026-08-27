"""Prepare and audit GPUMD qNEP extended-XYZ datasets with BEC labels."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Iterable

import numpy as np

from .bec_database import read_born, read_zborn


_PROPERTIES_RE = re.compile(r"(?i)(\bproperties\s*=\s*)(\"[^\"]*\"|\S+)")
_KEY_RE_TEMPLATE = r"(?i)(?:^|\s){key}\s*="
QNEP_BEC_OUTPUT_PRECISION = 10


@dataclass(frozen=True)
class ExtxyzFrame:
    natoms: int
    header: str
    atoms: tuple[str, ...]


@dataclass(frozen=True)
class BecData:
    tensors: np.ndarray
    labels: tuple[str, ...] = ()
    convention: str = "zstar"


def read_extxyz(path: str | Path) -> list[ExtxyzFrame]:
    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
    frames: list[ExtxyzFrame] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            natoms = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid atom count at {source}:{cursor + 1}") from exc
        if natoms < 1 or cursor + natoms + 1 >= len(lines):
            raise ValueError(f"Truncated extxyz frame at {source}:{cursor + 1}")
        frames.append(
            ExtxyzFrame(
                natoms=natoms,
                header=lines[cursor + 1],
                atoms=tuple(lines[cursor + 2:cursor + 2 + natoms]),
            )
        )
        cursor += natoms + 2
    if not frames:
        raise ValueError(f"No extxyz frames found in {source}")
    return frames


def _properties(header: str) -> tuple[re.Match[str], list[tuple[str, str, int]]]:
    match = _PROPERTIES_RE.search(header)
    if match is None:
        raise ValueError("Missing mandatory Properties field in extxyz header")
    value = match.group(2)
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    fields = value.split(":")
    if len(fields) % 3:
        raise ValueError(f"Malformed extxyz Properties field: {value}")
    properties: list[tuple[str, str, int]] = []
    for index in range(0, len(fields), 3):
        try:
            width = int(fields[index + 2])
        except ValueError as exc:
            raise ValueError(f"Invalid Properties width: {fields[index + 2]}") from exc
        properties.append((fields[index], fields[index + 1], width))
    return match, properties


def _property_offset(properties: Iterable[tuple[str, str, int]], name: str) -> tuple[int, int]:
    offset = 0
    for current, _kind, width in properties:
        if current.lower() == name.lower():
            return offset, width
        offset += width
    raise ValueError(f"Missing {name} property")


def _frame_labels(frame: ExtxyzFrame, properties: list[tuple[str, str, int]]) -> list[str]:
    offset, width = _property_offset(properties, "species")
    if width != 1:
        raise ValueError("species property must have width 1")
    expected = sum(item[2] for item in properties)
    labels: list[str] = []
    for atom_index, line in enumerate(frame.atoms, start=1):
        fields = line.split()
        if len(fields) != expected:
            raise ValueError(
                f"Frame atom {atom_index} has {len(fields)} columns; Properties declares {expected}"
            )
        labels.append(fields[offset])
    return labels


def _require_header_key(header: str, key: str) -> None:
    if not re.search(_KEY_RE_TEMPLATE.format(key=re.escape(key)), header):
        raise ValueError(f"Missing mandatory {key} field in extxyz header")


def read_bec_data(path: str | Path) -> BecData:
    """Read a canonical ZStar BEC source.

    Canonical tensors have displacement/force as rows and polarization/electric
    field as columns. qNEP conversion is deliberately done at export time.
    """

    source = Path(path)
    if source.suffix.lower() == ".json":
        data = json.loads(source.read_text(encoding="utf-8"))
        atoms = data.get("atoms")
        if not isinstance(atoms, list) or not atoms:
            raise ValueError(f"No atoms with BEC tensors in {source}")
        tensors = np.asarray([atom["tensor"] for atom in atoms], dtype=float)
        labels = tuple(str(atom.get("label", "")) for atom in atoms)
        convention = str(data.get("tensor_convention", "zstar"))
    elif source.name.upper().startswith("BORN"):
        _epsilon, tensors = read_born(source)
        labels = ()
        convention = "zstar"
    else:
        tensors = read_zborn(source)
        labels_list: list[str] = []
        for raw in source.read_text(encoding="utf-8", errors="ignore").splitlines():
            fields = raw.split()
            numeric = 0
            for field in fields:
                try:
                    float(field.replace("D", "E").replace("d", "e"))
                    numeric += 1
                except ValueError:
                    continue
            label_match = re.match(r"^\s*\*?\s*\d+\s+(\S+)", raw)
            if numeric >= 10 and label_match:
                labels_list.append(label_match.group(1))
        labels = tuple(labels_list) if len(labels_list) == len(tensors) else ()
        convention = "zstar"
    if tensors.ndim != 3 or tensors.shape[1:] != (3, 3):
        raise ValueError(f"BEC tensor array must have shape (natoms, 3, 3), got {tensors.shape}")
    return BecData(tensors=tensors, labels=labels, convention=convention)


def load_bec_map(path: str | Path) -> dict[int, Path]:
    source = Path(path).resolve()
    mapping: dict[int, Path] = {}
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = {"frame", "bec"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"BEC map missing columns: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            frame = int(row["frame"])
            if frame < 0:
                raise ValueError(f"Negative frame index at {source}:{row_number}")
            if frame in mapping:
                raise ValueError(f"Duplicate frame {frame} in {source}")
            bec = Path(row["bec"].strip())
            if not bec.is_absolute():
                bec = (source.parent / bec).resolve()
            mapping[frame] = bec
    return mapping


def _qnep_tensors(data: BecData) -> np.ndarray:
    convention = data.convention.lower()
    if "rows=atomic displacement" in convention or convention == "zstar":
        return np.transpose(data.tensors, (0, 2, 1))
    if "rows=electric" in convention or convention == "qnep":
        return data.tensors.copy()
    raise ValueError(f"Unknown BEC tensor convention: {data.convention}")


def augment_qnep_dataset(
    input_xyz: str | Path,
    output_xyz: str | Path,
    *,
    bec: str | Path | None = None,
    frame: int = 0,
    bec_map: str | Path | None = None,
    audit_output: str | Path | None = None,
) -> dict:
    """Append optional qNEP ``bec:R:9`` labels to selected extxyz frames."""

    if (bec is None) == (bec_map is None):
        raise ValueError("Provide exactly one of bec or bec_map")
    frames = read_extxyz(input_xyz)
    mapping = load_bec_map(bec_map) if bec_map is not None else {int(frame): Path(bec).resolve()}
    unknown = sorted(set(mapping) - set(range(len(frames))))
    if unknown:
        raise ValueError(f"BEC map refers to absent frame indices: {unknown}")

    output_lines: list[str] = []
    frame_audit: list[dict] = []
    for frame_index, current in enumerate(frames):
        _require_header_key(current.header, "lattice")
        _require_header_key(current.header, "energy")
        match, properties = _properties(current.header)
        _property_offset(properties, "pos")
        force_name = next((name for name, _kind, _width in properties if name.lower() in {"force", "forces"}), None)
        if force_name is None:
            raise ValueError(f"Frame {frame_index} has no force/forces property")
        labels = _frame_labels(current, properties)
        if any(name.lower() == "bec" for name, _kind, _width in properties):
            raise ValueError(f"Frame {frame_index} already contains a bec property")

        output_header = current.header
        output_atoms = list(current.atoms)
        entry = {
            "frame": frame_index,
            "natoms": current.natoms,
            "bec_labeled": False,
            "bec_source": None,
        }
        if frame_index in mapping:
            data = read_bec_data(mapping[frame_index])
            if len(data.tensors) != current.natoms:
                raise ValueError(
                    f"Frame {frame_index} has {current.natoms} atoms but {mapping[frame_index]} "
                    f"has {len(data.tensors)} BEC tensors"
                )
            if data.labels and tuple(labels) != data.labels:
                raise ValueError(
                    f"Atom order mismatch for frame {frame_index}: extxyz={labels}, "
                    f"BEC={list(data.labels)}"
                )
            qnep = _qnep_tensors(data)
            property_text = match.group(2)
            quoted = property_text.startswith('"')
            raw_properties = property_text[1:-1] if quoted else property_text
            new_properties = raw_properties + ":bec:R:9"
            if quoted:
                new_properties = f'"{new_properties}"'
            output_header = output_header[:match.start(2)] + new_properties + output_header[match.end(2):]
            output_atoms = [
                f"{line} "
                + " ".join(
                    f"{value:.{QNEP_BEC_OUTPUT_PRECISION}f}"
                    for value in tensor.reshape(9)
                )
                for line, tensor in zip(current.atoms, qnep)
            ]
            acoustic = np.sum(data.tensors, axis=0)
            entry.update(
                bec_labeled=True,
                bec_source=str(mapping[frame_index]),
                input_convention=data.convention,
                qnep_transform="transpose_zstar_to_field_rows",
                acoustic_sum_max_abs_e=float(np.max(np.abs(acoustic))),
            )
        output_lines.extend([str(current.natoms), output_header, *output_atoms])
        frame_audit.append(entry)

    target = Path(output_xyz).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(output_lines) + "\n", encoding="utf-8", newline="\n")
    summary = {
        "schema_version": 1,
        "format": "gpumd-qnep-extxyz",
        "input": str(Path(input_xyz).resolve()),
        "output": str(target),
        "frames": len(frames),
        "labeled_frames": sum(item["bec_labeled"] for item in frame_audit),
        "unlabeled_frames": sum(not item["bec_labeled"] for item in frame_audit),
        "bec_unit": "elementary_charge_e",
        "qnep_tensor_convention": "row-major rows=electric field/polarization; columns=force/displacement",
        "frames_detail": frame_audit,
    }
    audit = Path(audit_output).resolve() if audit_output else target.with_suffix(target.suffix + ".audit.json")
    audit.write_text(json.dumps(summary, indent=2), encoding="utf-8", newline="\n")
    summary["audit_output"] = str(audit)
    return summary


def check_qnep_dataset(path: str | Path, *, audit_output: str | Path | None = None) -> dict:
    frames = read_extxyz(path)
    details: list[dict] = []
    elements: set[str] = set()
    for index, frame in enumerate(frames):
        _require_header_key(frame.header, "lattice")
        _require_header_key(frame.header, "energy")
        _match, properties = _properties(frame.header)
        labels = _frame_labels(frame, properties)
        elements.update(labels)
        _property_offset(properties, "pos")
        if not any(name.lower() in {"force", "forces"} for name, _kind, _width in properties):
            raise ValueError(f"Frame {index} has no force/forces property")
        bec = next((item for item in properties if item[0].lower() == "bec"), None)
        if bec is not None and (bec[1].upper(), bec[2]) != ("R", 9):
            raise ValueError(f"Frame {index} BEC property must be bec:R:9")
        details.append({"frame": index, "natoms": frame.natoms, "has_bec": bec is not None})
    summary = {
        "schema_version": 1,
        "format": "gpumd-qnep-extxyz",
        "path": str(Path(path).resolve()),
        "frames": len(frames),
        "labeled_frames": sum(item["has_bec"] for item in details),
        "elements": sorted(elements),
        "valid": True,
        "frames_detail": details,
    }
    if audit_output:
        target = Path(audit_output).resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(summary, indent=2), encoding="utf-8", newline="\n")
        summary["audit_output"] = str(target)
    return summary


def write_qnep_input(
    dataset: str | Path,
    output: str | Path = "nep.in",
    *,
    charge_mode: int = 2,
    lambda_z: float = 0.5,
) -> Path:
    if charge_mode not in {1, 2}:
        raise ValueError("qNEP charge_mode must be 1 or 2")
    if lambda_z < 0:
        raise ValueError("lambda_z must be nonnegative")
    summary = check_qnep_dataset(dataset)
    if summary["labeled_frames"] == 0:
        raise ValueError("Dataset has no bec:R:9 labels")
    elements = summary["elements"]
    target = Path(output).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        "# Minimal qNEP input generated by ZStar; converge model capacity and training controls.\n"
        f"type {len(elements)} {' '.join(elements)}\n"
        f"charge_mode {charge_mode}\n"
        f"lambda_z {lambda_z:.10g}\n",
        encoding="utf-8",
        newline="\n",
    )
    return target
