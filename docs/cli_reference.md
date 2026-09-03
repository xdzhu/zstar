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
| `zstar bec` | `pre/run/job/stat/post` | Polarization, APT/BEC, and `BORN`; calculators: ABACUS + PYATB, VASP, CP2K, QE. |
| `zstar phonon` (`ph`) | `pre/run/job/stat/post/irrep` | Displacements, serial force calculations, force constants, frequencies, and Gamma irreps. |
| `zstar spectra` | `pre/run/job/stat/post` | IR and Raman workflows for ABACUS + PYATB, VASP, CP2K, and QE. |
| `zstar dielectric` (`diel`) | `static` (`zero`), `freq`, `optics` | Ionic static response, frequency-dependent vibrational response, and electronic optics. |
| `zstar backend list` | `--check`, `--json`, `--discover` | List implemented capabilities and optionally check configured executables or plugins. |
| `zstar config` | `init/show/set/check` | Layered executable and launch configuration. |
| `zstar response` | `validate/import-bec/import-abacus/import-phonopy/intrinsic` | Calculator-neutral response documents and intrinsic low-dimensional response. |
| `zstar density` | `vasp-cube/qe-input/qe-sidecar/cp2k-block/sidecar` | Density-export adapters and provenance sidecars. |
| `zstar stru` | `convert/wyckoff` | Structure conversion and symmetry inspection. |
| `zstar data` | `qnep/db` | qNEP training-data export and traceable BEC/High-K databases. |
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
Environment setup belongs in the generated job's `--env-script`.

## Representative lifecycles

Run the following examples from a prepared workflow directory containing the
calculator inputs. ABACUS + PYATB bulk BEC:

```bash
zstar bec pre --stru STRU
zstar bec job --root . --system slurm --tasks 28 --env-script env.sh
zstar bec stat --root .
zstar bec post --root .
```

For `zstar bec pre`, ABACUS + PYATB and the forward finite-difference method
are the defaults. The calculator and `--pyatb` switch are therefore omitted
from the canonical example; specify `--calculator cp2k`, `vasp`, or `qe` only
when changing backends.

Phonons and mode classification:

```bash
zstar phonon pre --root . --calculator abacus --stru STRU --dim "2 2 2"
zstar phonon job --root . --system slurm --tasks 28
zstar phonon stat --root .
zstar phonon post --root . --nac
zstar phonon irrep --root . --file irreps.yaml --mode db
```

IR and Raman:

```bash
zstar spectra pre --calculator abacus --kind all --root spectra \
  --stru STRU --qpoints qpoints.yaml --born Z-BORN-symm.out --dielectric BORN
zstar spectra job --root spectra --system slurm --tasks 28
zstar spectra stat --root spectra
zstar spectra post --root spectra
```

Static and frequency-dependent dielectric response:

```bash
zstar dielectric static --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 3
zstar dielectric freq --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 3
```

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
| `zstar spectra pre/run/job/stat/post` | `prepare/run/status/collect/script`; low-level `ir` and `raman` remain available |
| `zstar dielectric static/freq` | `zstar calc/freq` |
| `zstar stru convert/wyckoff` | `zstar vasp/wyckoff` |
| `zstar skill install/path/preflight` | `zstar agent-skill ...` |
