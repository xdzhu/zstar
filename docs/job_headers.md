# Job headers and calculation environments

[Chinese tutorial](job_headers.zh-CN.md)

ZStar generates one serial, resumable driver for a calculation directory.
The job header controls the scheduler allocation and environment setup. The
configuration controls executable paths and MPI ranks/OMP threads. Generating
a script does not submit it.

## Header selection

Exactly one header is embedded, in this order:

| Priority | Name | Location |
| --- | --- | --- |
| 1 | Specified | `--header /path/to/header.sh` |
| 2 | Current | `header.sh` in the workflow root (`--root`, default `.`) |
| 3 | Global | `~/.zstar/header.sh` |
| 4 | Default | Generated, editable scheduler template |

Headers are not merged, and parent directories are not searched. An explicitly
specified missing file is an error, not a request to fall back. A selected empty
file is also an error. Use a comment-only file for an intentionally empty header.
The selected header replaces automatic scheduler directives, including legacy
`--queue`, `--account`, `--nodes` and `--walltime` defaults. Record executable
paths and parallelism separately:

```bash
zstar config set executables.abacus /opt/abacus/bin/abacus --user
zstar config set execution.mpi 1 --user
zstar config set execution.omp 40 --user
zstar config set abacus.pseudo_dir /data/PSEUDO --user
zstar config set abacus.orbital_dir /data/ORBITAL --user
```

Omit `--user` for workflow-local configuration. Existing `execution.tasks` and
`execution.cpus_per_task` keys remain supported; when both forms are present,
`mpi` and `omp` take precedence. Explicit `--tasks` / `--cpus-per-task` overrides
remain supported. MPI/OMP must fit the header allocation; ZStar does not infer
launch settings by parsing arbitrary shell code or scheduler directives.

## Slurm example

Put the following in `header.sh` in the workflow root (or in
`~/.zstar/header.sh` to reuse it globally):

```bash
#!/usr/bin/env bash
#SBATCH --job-name=zstar
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --output=zstar-%j.out

module load compiler mpi
source /opt/conda/etc/profile.d/conda.sh
conda activate zstar
```

Replace site-specific names and paths before use. All scheduler directives must
precede the first shell command. Generate and submit from the workflow directory:

```bash
zstar bec pre --stru STRU
zstar bec job --system slurm
sbatch run_zstar_born.slurm
```

The same `--header` selection is available for `zstar phonon job` and
`zstar spectra job`, and for ABACUS, CP2K, VASP and QE job generators. For a
different cluster, choose a different header explicitly:

```bash
zstar bec job --system slurm --header /data/cluster-b/header.sh
```

## Torque and local shell

A Torque header for the same 1 MPI x 40 OMP allocation:

```bash
#!/usr/bin/env bash
#PBS -N zstar
#PBS -q compute
#PBS -l nodes=1:ppn=40
#PBS -l walltime=24:00:00

module load compiler mpi
source /opt/conda/etc/profile.d/conda.sh
conda activate zstar
```

```bash
zstar bec job --system torque
qsub run_zstar_born.pbs
```

For `--system shell`, the header contains only environment commands, without
`#PBS` or `#SBATCH` directives; run `bash run_zstar_born.sh`. A header for a
different scheduler is rejected. This prevents a global Slurm header from
silently becoming a local-shell environment script. Pass `--header` to switch.

## Reproducibility and restart

The generated script embeds the selected header, rather than sourcing that file
when the job starts. `.zstar/job_header.json` records its path, SHA-256 and
selection level. Editing the original header does not change an already
generated script; regenerate it when the environment changes. Scripts use Bash
strict error handling and run header commands from the workflow root.

Check the selected script before submission, especially allocation, account,
module names and MPI/OMP. The stage-state files remain in `.zstar`; rerunning
the same driver resumes unfinished stages. `zstar bec stat` reports progress.
Header support does not change reference-first execution or restart-density
copying. Legacy `--env-script` remains supported for existing workflows; place
new module/source commands directly in a header.

## No header configured

ZStar still produces a complete default script. Its compact comment block shows
Specified, Current and Global locations and links to this tutorial. Edit its
scheduler lines directly for a one-off job, or create a reusable header before
generating the next script. Phonopy and PYATB must be installed in the activated
Python environment; DFT executables can be found on `PATH` or configured above.
