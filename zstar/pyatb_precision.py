"""Run PYATB with a lossless polarization writer, without editing PYATB files.

The numerical kernels and their constants are unchanged. This narrowly scoped
adapter preserves the original six-decimal file alongside the full-precision
output needed for small mixed-direction finite displacements.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex
import shutil
import sys

import numpy as np


def write_precise_polarization(response):
    target = Path(response.output_path) / "polarization.dat"
    original = target.read_bytes()
    arrays = [np.asarray(getattr(response, name), dtype=float) for name in
              ("polarization_ion", "polarization_ele", "polarization", "modulus")]
    if any(a.shape != (3,) or not np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("PYATB polarization writer API changed or produced invalid values")
    ion, ele, polar, quantum = arrays
    lines = ["The Ionic Phase      : " + " ".join(f"{v:.16e}" for v in ion),
             "The Electronic Phase : " + " ".join(f"{v:.16e}" for v in ele)]
    lines += [f"The calculated polarization direction is in {axis}, P = {p:.16e} (mod {q:.16e}) C/m^2."
              for axis, p, q in zip("abc", polar, quantum)]
    target.with_name("polarization.rounded.dat").write_bytes(original)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    target.with_name("zstar_precision.json").write_text(json.dumps({
        "adapter": "zstar.pyatb_precision", "numerical_kernel_changed": False,
        "original_sha256": hashlib.sha256(original).hexdigest(),
        "precision_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "format": "16 digits after decimal point in scientific notation",
    }, indent=2) + "\n", encoding="utf-8")


def precision_command(command, executable="pyatb"):
    if "-m zstar.pyatb_precision" in command:
        return command
    tokens = shlex.split(command)
    if any(t in (";", "&&", "||", "|", ">", "<") for t in tokens):
        raise ValueError("For a shared ensemble use a direct PYATB launcher command, or launch python -m zstar.pyatb_precision inside your wrapper")
    matches = [i for i, t in enumerate(tokens) if t == executable or Path(t).name == "pyatb"]
    if len(matches) != 1:
        raise ValueError("Shared mixed displacements require full-precision polarization. Set the PYATB command to 'mpirun -np N python -m zstar.pyatb_precision' in the PYATB/ZStar environment.")
    index = matches[0]
    binary = shutil.which(tokens[index])
    python = sys.executable
    if binary:
        with Path(binary).open("rb") as handle:
            first = handle.readline(4096).decode("utf-8", errors="replace").strip()
        if first.startswith("#!"):
            interpreter = shlex.split(first[2:])
            if len(interpreter) == 1 and Path(interpreter[0]).is_file():
                python = interpreter[0]
    tokens[index:index + 1] = [python, "-m", "zstar.pyatb_precision"]
    return shlex.join(tokens)


def main():
    import pyatb
    # New PYATB initializes lazily; imported response modules capture its
    # runtime globals. Older releases initialize at package import instead.
    initialize = getattr(pyatb, "initialize_runtime", None)
    if initialize is not None:
        initialize()
    from pyatb.berry.polarization import Polarization
    from pyatb.main import main as pyatb_main

    original = getattr(Polarization, "print_data", None)
    if not callable(original):
        raise RuntimeError("Unsupported PYATB polarization writer API; no numerical calculation was changed")

    def precise(self):
        original(self)
        write_precise_polarization(self)

    Polarization.print_data = precise
    try:
        return pyatb_main()
    finally:
        Polarization.print_data = original


if __name__ == "__main__":
    sys.exit(main())
