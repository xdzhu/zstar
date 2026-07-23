# Scheduler Drivers

The recommended Born-charge driver is generated for the active environment
instead of being copied into every displacement directory:

```bash
# Local or direct compute node
zstar workflow script --backend shell \
  --tasks 32 --cpus-per-task 1 --env-script /path/to/zstar-env.sh

# Slurm
zstar workflow script --backend slurm \
  --nodes 1 --tasks 32 --cpus-per-task 1 \
  --queue compute --walltime 24:00:00 \
  --env-script /path/to/zstar-env.sh

# Torque/PBS
zstar workflow script --backend torque \
  --nodes 1 --tasks 32 --cpus-per-task 1 \
  --queue batch --walltime 24:00:00 \
  --env-script /path/to/zstar-env.sh
```

The generated script runs `0.no-move` first, applies the reference insulation
gate, reuses the converged reference charge density, executes displacement
stages serially, and resumes from `.zstar/stages/`.

Backend-aware launch defaults are:

| Backend | ABACUS/PyATB launcher | Generated file |
| --- | --- | --- |
| Shell | `mpirun -np N` | `run_zstar_born.sh` |
| Slurm | `srun --ntasks=N` | `run_zstar_born.slurm` |
| Torque/PBS | `mpirun -np N` | `run_zstar_born.pbs` |

Use `--abacus-command` or `--pyatb-command` to override a launcher. Before a
production run, generate with `--dry-run` to check the environment, ordering,
state files, and scheduler directives without invoking either solver.

`job_BORN.slurm` is retained as a legacy site-specific example. New workflows
should use `zstar workflow script`, which is portable across the three
supported backends and records its choices in
`.zstar/backend_manifest.json`.
