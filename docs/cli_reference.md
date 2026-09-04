# ZStar command-line reference

ZStar groups public commands by scientific object. Calculation workflows use
the same lifecycle wherever the operation exists:

```text
pre -> run -> job -> stat -> post
```

- `pre` creates inputs and writes `.zstar/<family>.json`.
- `run` executes the prepared stages serially and resumes completed stages.
- `job` writes one shell, Slurm, or Torque driver for the complete serial chain.
- `stat` reports prepared, running, completed, failed, and blocked stages.
- `post` validates and collects scientific outputs.

`prepare`, `status`, `collect`/`deal`, and `script` remain accepted compatibility
aliases. New documentation and automation should use the short canonical verbs.

## Public command tree

| Family | Canonical actions | Purpose |
| --- | --- | --- |
| `zstar bec` | `pre/job/run/stat/post` | Polarization, APT/BEC, and `BORN`; calculators: ABACUS + PYATB, VASP, CP2K, QE. |
| `zstar phonon` (`ph`) | `pre/job/run/stat/post/irrep` | Displacements, serial force calculations, force constants, frequencies, and Gamma irreps. |
| `zstar spectra` | `pre/job/run/stat/post` | IR and Raman workflows for ABACUS + PYATB, VASP, CP2K, and QE. |
| `zstar dielectric` (`diel`) | `static` (`zero`), `freq`, `optics` | Ionic static response, frequency-dependent vibrational response, and electronic optics. |
| `zstar backend list` | `--check`, `--json`, `--discover` | List implemented capabilities and optionally check configured executables or plugins. |
| `zstar config` | `init/show/set/check` | Layered executable and launch configuration. |
| `zstar response` | `validate/import-bec/import-abacus/import-phonopy/intrinsic` | Calculator-neutral response documents and intrinsic low-dimensional response. |
| `zstar density` | `vasp-cube/qe-input/qe-sidecar/cp2k-block/sidecar` | Density-export adapters and provenance sidecars. |
| `zstar stru` | `convert/wyckoff` | Structure conversion and symmetry inspection. |
| `zstar data` | `db/qnep` | Traceable BEC/High-K databases and qNEP training-data export. |
| `zstar skill` | `install/path/preflight` | Install or inspect the packaged Agent Skill and run non-mutating preflight checks. |
| `zstar pot` | option-driven | Axis profiles, plane maps, directional profiles, vacuum steps, and mirror asymmetry. |

The old fine-grained commands (`gen`, `workflow`, `deal`, `postph`, `ir`,
`raman`, calculator-specific `*-bec`, and others) remain available as a
compatibility and expert layer. `zstar backend list` is the only public backend
entry.

## Calculator configuration

Initialize a project-local configuration and edit executable paths once:

```bash
zstar config init
zstar config set executables.abacus /opt/abacus/bin/abacus
zstar config set executables.pyatb /opt/pyatb/bin/pyatb
zstar config set executables.vasp /opt/vasp/bin/vasp_std
zstar config set executables.cp2k /opt/cp2k/bin/cp2k.psmp
zstar config check
```

The supported keys are `abacus`, `pyatb`, `vasp`, `cp2k`, `qe_pw`, `qe_ph`,
`qe_dynmat`, and `phonopy`. Resolution order, from lowest to highest priority,
is built-in defaults, user config, project config, and environment variables.
The project file is `.zstar/config.toml`; the user file is
`~/.config/zstar/config.toml` on Linux and `%APPDATA%/zstar/config.toml` on
Windows. `ZSTAR_CONFIG` selects an alternative user file. A one-run command-line
override has the highest priority.

Environment variables use names such as `ZSTAR_ABACUS_EXECUTABLE` and
`ZSTAR_QE_PW_EXECUTABLE`. Store the calculator executable in configuration;
ZStar adds `mpirun -np N` for shell/Torque jobs or `srun --ntasks=N` for Slurm.
Environment setup belongs in the selected job header.

### Cluster resources and environment

Queue, account, resources, time limits and module/source commands belong in a
single selected header: **Specified** `--header FILE`, then **Current**
`header.sh` in the workflow root, then **Global** `~/.zstar/header.sh`. There is
no merge or parent-directory search. If none exists, an editable default is
generated. MPI/OMP and executable paths remain configuration settings:

```bash
zstar config set execution.mpi 1
zstar config set execution.omp 40
zstar bec job --system slurm --header /path/to/header.sh
```

The same selection applies to `zstar phonon job` and `zstar spectra job`.
Selected contents are embedded in the script and hashed in
`.zstar/job_header.json`. Keep MPI/OMP consistent with the scheduler allocation.
Legacy resource flags and `--env-script` remain supported. For complete header
examples, default behavior, execution and restart, see [the tutorial](job_headers.md).

### ABACUS pseudopotentials and orbitals

`zstar bec pre` resolves ABACUS resources before creating displacement folders.
Use explicit directories for a case that differs from the user's defaults:

```bash
zstar bec pre --stru STRU \
  --pp /path/to/PSEUDO \
  --orb /path/to/ORBITAL
```

Otherwise configure the global defaults with:

```bash
zstar config set abacus.pseudo_dir /data/PSEUDO --user
zstar config set abacus.orbital_dir /data/ORBITAL --user
```

ZStar first tries the exact filename from `STRU`, then accepts a unique
element-prefix match such as `Si_*.upf` or `Si_*.orb`.  Ambiguous matches stop
with a candidate list and explain how to select an exact file or narrower
directory.  The input `STRU` is left unchanged; `.zstar/STRU.resolved` and
`.zstar/assets.json` provide the generated copy and resource provenance.

## Representative lifecycles

Run the following examples from a prepared workflow directory containing the
calculator inputs. ABACUS + PYATB bulk BEC:

```bash
zstar bec pre --stru STRU
zstar bec job --system slurm --header header.sh
zstar bec stat --root .
zstar bec post --root .
```

In `0.3.0rc2`, `zstar bec pre` defaults to ABACUS + PYATB and
Phonopy's symmetry-adapted Unified BEC/Gamma displacement ensemble, with
automatic +/- selection. The calculator and `--pyatb` switch remain optional;
specify `--calculator cp2k`, `vasp`, or `qe` only when changing backends.
`--ensemble cartesian` retains the legacy atom/direction layout. See the
[shared response tutorial](research/shared_response/USAGE.md) for scope,
actual displacement units, raw diagnostics, and compatibility.

For the shared Gamma route, `zstar bec post` already generates the phonon
outputs. Prepare finite-q/supercell phonons in a separate directory:

```bash
zstar phonon pre --stru STRU --dim "2 2 2"
zstar phonon job --system slurm --tasks 28
zstar phonon stat
zstar phonon post
zstar phonon irrep
```

IR and Raman:

```bash
zstar spectra pre --stru STRU
zstar spectra job --system slurm --tasks 28
zstar spectra stat
zstar spectra post
```

Static and frequency-dependent dielectric response:

```bash
zstar dielectric static
zstar dielectric freq
```

Low-dimensional outputs use their documented source-field normalization.
For independently screened macroscopic slab data only, the optional
`--slab-boundary macroscopic` with `--thickness` converts the total tensor
using direct in-plane and inverse out-of-plane response. It does not add
missing local-field physics to PYATB results. See
[response conventions](response_conventions.md) and the
[eight-system benchmarks](../examples/Shared_Response/README.md).

## Electrostatic-potential coverage

`zstar pot` retains the complete analysis used in the ZStar paper:

```bash
zstar pot --cube ElecStaticPot.cube --axes z \
  --plane xy --plane-average --tile 5 5 \
  --direction a+b --mirror-test \
  --vacuum-sides --vacuum-window 0.75 --outdir potential
```

This command writes the one-dimensional potential, a tiled plane map with a
dashed central unit cell, a perpendicular-plane-averaged directional profile,
optimized one-period mirror analysis, and the two-surface vacuum potential
step. The mirror metric tests asymmetry within one selected periodic direction;
it does not compare `a+b` against `a-b`.

## Compatibility map

| Canonical command | Accepted legacy form |
| --- | --- |
| `zstar bec pre` | `zstar gen` |
| `zstar bec run/job/stat` | `zstar workflow run/status/script` |
| `zstar bec post` | `zstar deal` |
| `zstar phonon pre/post/irrep` | `zstar ph/postph/irrep` |
| `zstar spectra pre/job/run/stat/post` | `prepare/script/run/status/collect`; low-level `ir` and `raman` remain available |
| `zstar dielectric static/freq` | `zstar calc/freq` |
| `zstar stru convert/wyckoff` | `zstar vasp/wyckoff` |
| `zstar skill install/path/preflight` | `zstar agent-skill ...` |
