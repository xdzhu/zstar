"""Compatibility helpers for old and new PYATB optical-response interfaces."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Optional, Sequence, Tuple

import numpy as np


DEFAULT_LEGACY_OMEGA_MAX_EV = 30.0
DEFAULT_LEGACY_DOMEGA_EV = 0.1


_FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")
_OPTICAL_BLOCK_RE = re.compile(
    r"(?ms)^[ \t]*OPTICAL_CONDUCTIVITY[ \t]*\n[ \t]*\{.*?^[ \t]*\}[ \t]*$"
)


@dataclass(frozen=True)
class PyATBCapabilities:
    """Features exposed by the PYATB installation used for a workflow."""

    version: str
    executable: str
    static_dielectric_only: bool
    detection: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class BandGapResult:
    gap_eV: float
    vbm_eV: Optional[float]
    cbm_eV: Optional[float]
    insulating: bool
    threshold_eV: float
    source: str

    def to_dict(self) -> dict:
        return asdict(self)


def _current_installation_has_static() -> bool:
    spec = importlib.util.find_spec("pyatb")
    if spec is None or spec.origin is None:
        return False
    root = Path(spec.origin).resolve().parent
    candidates = (
        root / "io" / "default_input.py",
        root / "berry" / "optical_conductivity.py",
    )
    for path in candidates:
        try:
            if "static_dielectric_only" in path.read_text(
                encoding="utf-8", errors="ignore"
            ):
                return True
        except OSError:
            continue
    return False


def _external_installation_probe(executable: str) -> Optional[Tuple[str, bool]]:
    """Probe the Python environment named by a console-script shebang."""

    path = Path(executable)
    try:
        first_line = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    except (OSError, IndexError):
        return None
    if not first_line.startswith("#!"):
        return None
    python_executable = first_line[2:].strip().split()[0]
    if not Path(python_executable).exists():
        return None

    probe = (
        "import importlib.metadata as m, importlib.util, pathlib;"
        "s=importlib.util.find_spec('pyatb');"
        "r=pathlib.Path(s.origin).resolve().parent if s and s.origin else None;"
        "ps=[r/'io'/'default_input.py',r/'berry'/'optical_conductivity.py'] "
        "if r else [];"
        "print(m.version('pyatb'));"
        "print(int(any(p.exists() and 'static_dielectric_only' in "
        "p.read_text(errors='ignore') for p in ps)))"
    )
    try:
        result = subprocess.run(
            [python_executable, "-c", probe],
            check=True,
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    lines = result.stdout.strip().splitlines()
    if len(lines) < 2:
        return None
    return lines[-2].strip(), lines[-1].strip() == "1"


def _adjacent_target_installation_probe(
    executable: str,
) -> Optional[Tuple[str, bool]]:
    """Detect a ``pip --target DIR`` install with console scripts in ``DIR/bin``."""

    path = Path(executable).resolve()
    target = path.parent.parent
    package = target / "pyatb"
    if path.parent.name != "bin" or not package.is_dir():
        return None
    candidates = (
        package / "io" / "default_input.py",
        package / "berry" / "optical_conductivity.py",
    )
    has_static = any(
        source.is_file()
        and "static_dielectric_only"
        in source.read_text(encoding="utf-8", errors="ignore")
        for source in candidates
    )
    versions = sorted(target.glob("pyatb-*.dist-info"))
    version = (
        versions[-1].name[len("pyatb-") : -len(".dist-info")]
        if versions
        else "unknown"
    )
    return version, has_static


def detect_pyatb_capabilities(executable: str = "pyatb") -> PyATBCapabilities:
    """Detect whether the selected PYATB accepts ``static_dielectric_only``."""

    resolved = shutil.which(executable) or executable
    adjacent = _adjacent_target_installation_probe(resolved)
    if adjacent is not None:
        version, has_static = adjacent
        return PyATBCapabilities(
            version=version,
            executable=resolved,
            static_dielectric_only=has_static,
            detection="adjacent pip-target installation",
        )
    external = _external_installation_probe(resolved)
    if external is not None:
        version, has_static = external
        return PyATBCapabilities(
            version=version,
            executable=resolved,
            static_dielectric_only=has_static,
            detection="console-script environment",
        )

    try:
        version = importlib.metadata.version("pyatb")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return PyATBCapabilities(
        version=version,
        executable=resolved,
        static_dielectric_only=_current_installation_has_static(),
        detection="current Python environment",
    )


def _replace_or_append_parameter(block: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^([ \t]*){re.escape(key)}\s+.*$")
    match = pattern.search(block)
    if match:
        indent = match.group(1)
        return pattern.sub(f"{indent}{key:<28}{value}", block, count=1)
    close = block.rfind("}")
    if close < 0:
        raise ValueError("Malformed OPTICAL_CONDUCTIVITY block: missing closing brace")
    return block[:close] + f"    {key:<28}{value}\n" + block[close:]


def _extract_scalar(text: str, key: str) -> Optional[int]:
    match = re.search(rf"(?m)^[ \t]*{re.escape(key)}[ \t]+(\d+)", text)
    return int(match.group(1)) if match else None


def configure_optical_input(
    input_path: str | Path,
    *,
    capabilities: Optional[PyATBCapabilities] = None,
    static_only: bool = True,
    legacy_omega_max: float = DEFAULT_LEGACY_OMEGA_MAX_EV,
    legacy_domega: float = DEFAULT_LEGACY_DOMEGA_EV,
    eta: float = 0.05,
    grid: Optional[Sequence[int]] = None,
) -> dict:
    """Configure a compact, converged electronic dielectric calculation.

    New PYATB builds use their direct static-response kernel. Older builds are
    evaluated on a coarse optical grid through ``legacy_omega_max``. The
    default 0--30 eV range with 0.1 eV spacing was selected by comparing the
    legacy spectrum intercept with the direct static response.
    """

    path = Path(input_path)
    text = path.read_text(encoding="utf-8")
    caps = capabilities or detect_pyatb_capabilities()
    use_direct_static = bool(static_only and caps.static_dielectric_only)

    match = _OPTICAL_BLOCK_RE.search(text)
    if match:
        block = match.group(0)
    else:
        occ_band = _extract_scalar(text, "occ_band")
        if occ_band is None:
            raise ValueError(
                "Cannot add OPTICAL_CONDUCTIVITY: no occ_band was found in the PYATB input"
            )
        if grid is None:
            nk = [_extract_scalar(text, f"nk{i}") for i in (1, 2, 3)]
            grid_values = [value if value is not None else 1 for value in nk]
        else:
            grid_values = [int(value) for value in grid]
        if len(grid_values) != 3 or any(value < 1 for value in grid_values):
            raise ValueError("Optical integration grid must contain three positive integers")
        block = (
            "OPTICAL_CONDUCTIVITY\n"
            "{\n"
            f"    {'occ_band':<28}{occ_band}\n"
            f"    {'omega':<28}0.0 {legacy_omega_max:g}\n"
            f"    {'domega':<28}{legacy_domega:g}\n"
            f"    {'eta':<28}{eta:g}\n"
            f"    {'grid':<28}{grid_values[0]} {grid_values[1]} {grid_values[2]}\n"
            "    method                      1\n"
            "}"
        )

    block = _replace_or_append_parameter(
        block, "omega", f"0.0 {float(legacy_omega_max):g}"
    )
    block = _replace_or_append_parameter(block, "domega", f"{float(legacy_domega):g}")
    block = _replace_or_append_parameter(block, "eta", f"{float(eta):g}")
    if use_direct_static:
        block = _replace_or_append_parameter(block, "static_dielectric_only", "1")
    else:
        block = re.sub(
            r"(?m)^[ \t]*static_dielectric_only\s+.*\n?", "", block
        )

    if match:
        output = text[: match.start()] + block + text[match.end() :]
    else:
        output = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(output, encoding="utf-8")

    report = {
        "input": str(path.resolve()),
        "pyatb": caps.to_dict(),
        "mode": "direct-static" if use_direct_static else "legacy-compact-spectrum",
        "legacy_omega_max_eV": float(legacy_omega_max),
        "legacy_domega_eV": float(legacy_domega),
    }
    report_path = path.parent / "zstar_pyatb_compat.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def _numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values = [float(value) for value in _FLOAT_RE.findall(stripped)]
        if values:
            rows.append(values)
    return rows


def read_static_dielectric(
    path_or_directory: str | Path,
) -> Tuple[np.ndarray, Path]:
    """Read a static 3x3 dielectric tensor from either PYATB output format."""

    root = Path(path_or_directory)
    if root.is_file():
        candidates = [root]
    else:
        candidates = [
            root / "static_dielectric_function.dat",
            root / "Out" / "Optical_Conductivity" / "static_dielectric_function.dat",
            root / "Optical_Conductivity" / "static_dielectric_function.dat",
            root / "dielectric_function_real_part.dat",
            root / "Out" / "Optical_Conductivity" / "dielectric_function_real_part.dat",
            root / "Optical_Conductivity" / "dielectric_function_real_part.dat",
        ]

    for candidate in candidates:
        if not candidate.is_file() or candidate.stat().st_size == 0:
            continue
        rows = _numeric_rows(candidate)
        if not rows:
            continue
        values = rows[0]
        if candidate.name == "dielectric_function_real_part.dat":
            if len(values) < 10:
                continue
            values = values[1:10]
        elif len(values) >= 9:
            values = values[-9:]
        else:
            continue
        tensor = np.asarray(values, dtype=float).reshape(3, 3)
        return tensor, candidate.resolve()
    raise FileNotFoundError(
        f"No PYATB static dielectric output found under {root.resolve()}"
    )


def read_band_gap(
    path_or_directory: str | Path,
    *,
    threshold_eV: float = 0.01,
) -> BandGapResult:
    """Read a PYATB band gap and apply an insulating threshold."""

    root = Path(path_or_directory)
    if root.is_file():
        candidates = [root]
    else:
        candidates = [
            root / "band_info.dat",
            *sorted(root.glob("Out/**/band_info.dat")),
        ]
    source = next(
        (
            path
            for path in candidates
            if path.is_file() and path.stat().st_size > 0
        ),
        None,
    )
    if source is None:
        raise FileNotFoundError(
            f"No PYATB band_info.dat found under {root.resolve()}"
        )
    text = source.read_text(encoding="utf-8", errors="ignore")

    def extract(label: str) -> Optional[float]:
        match = re.search(
            rf"{re.escape(label)}\s*({_FLOAT_RE.pattern})", text, re.IGNORECASE
        )
        return float(match.group(1)) if match else None

    gap = extract("Band gap (eV):")
    vbm = extract("Eigenvalue of VBM (eV):")
    cbm = extract("Eigenvalue of CBM (eV):")
    if gap is None:
        # PYATB omits VBM/CBM data when no insulating separation is found.
        gap = 0.0
    return BandGapResult(
        gap_eV=float(gap),
        vbm_eV=vbm,
        cbm_eV=cbm,
        insulating=bool(float(gap) >= float(threshold_eV)),
        threshold_eV=float(threshold_eV),
        source=str(source.resolve()),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Small diagnostic CLI used by generated workflow scripts."""

    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detect", action="store_true")
    parser.add_argument("--configure", metavar="INPUT")
    parser.add_argument("--read", metavar="OUTPUT")
    parser.add_argument("--pyatb", default="pyatb")
    parser.add_argument(
        "--legacy-omega-max",
        type=float,
        default=DEFAULT_LEGACY_OMEGA_MAX_EV,
    )
    args = parser.parse_args(argv)

    caps = detect_pyatb_capabilities(args.pyatb)
    if args.configure:
        report = configure_optical_input(
            args.configure,
            capabilities=caps,
            legacy_omega_max=args.legacy_omega_max,
        )
        print(json.dumps(report, indent=2))
    elif args.read:
        tensor, source = read_static_dielectric(args.read)
        print(source)
        np.savetxt(sys.stdout, tensor, fmt="%15.8e")
    else:
        print(json.dumps(caps.to_dict(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
