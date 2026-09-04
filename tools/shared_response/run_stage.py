"""Time one isolated ABACUS stage using a monotonic clock and explicit cores."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import socket
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", type=Path)
    parser.add_argument("--mpi", type=int, default=1)
    parser.add_argument("--omp", type=int, default=40)
    parser.add_argument("--timeout", type=float, default=14400)
    args = parser.parse_args()
    stage = args.stage.resolve()
    stamp = stage / "timing.json"
    if stamp.exists():
        raise SystemExit("timing.json already exists; use a new run directory")
    command = ["mpirun", "-np", str(args.mpi),
               "/home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus"]
    env = {**os.environ, "OMP_NUM_THREADS": str(args.omp),
           "MKL_NUM_THREADS": str(args.omp), "OPENBLAS_NUM_THREADS": "1"}
    record = {"stage": str(stage), "host": socket.gethostname(),
              "mpi": args.mpi, "omp": args.omp, "allocated_cores": args.mpi * args.omp,
              "command": command, "started_utc": datetime.now(timezone.utc).isoformat(),
              "pid": os.getpid(), "status": "running"}
    stamp.write_text(json.dumps(record, indent=2) + "\n")
    start = time.monotonic()
    try:
        with (stage / "abacus.stdout").open("w") as out, (stage / "abacus.stderr").open("w") as err:
            process = subprocess.run(command, cwd=stage, env=env, stdout=out,
                                     stderr=err, timeout=args.timeout)
        record["returncode"] = process.returncode
        record["status"] = "process_finished" if process.returncode == 0 else "failed"
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = repr(exc)
        raise
    finally:
        record["wall_seconds"] = time.monotonic() - start
        record["allocated_core_hours"] = record["wall_seconds"] * args.mpi * args.omp / 3600
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        stamp.write_text(json.dumps(record, indent=2) + "\n")
    raise SystemExit(record["returncode"])


if __name__ == "__main__":
    main()
