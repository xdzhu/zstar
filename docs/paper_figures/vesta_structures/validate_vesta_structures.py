"""Validate exported VESTA structures against their ABACUS sources with ASE."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.io import read


def main() -> None:
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "manifest.json").read_text())
    failures: list[str] = []
    for record in manifest["structures"]:
        atoms = read(root / record["output"], format="vasp")
        if len(atoms) != record["atoms"]:
            failures.append(f"{record['output']}: atom count")
        volume = abs(float(np.linalg.det(atoms.cell.array)))
        if not np.isclose(volume, record["cell_volume_A3"], rtol=1.0e-9, atol=1.0e-8):
            failures.append(f"{record['output']}: cell volume")
        scaled = atoms.get_scaled_positions(wrap=False)
        if not np.isfinite(scaled).all():
            failures.append(f"{record['output']}: non-finite coordinates")
        print(f"OK {record['output']}: {len(atoms)} atoms, {volume:.6f} A^3")
    if failures:
        raise SystemExit("Validation failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
