# Agent-Friendly ZStar Workflows

ZStar ships the standards-compliant Agent Skill `run-zstar-workflows`. The skill
turns the command-line interface, physical conventions, resumable state, and
output contracts into reusable instructions that a compatible coding or
scientific agent can discover and follow.

The skill does not replace scientific judgment or grant permission to launch
calculations. It helps an agent choose the correct workflow, inspect readiness,
retain provenance, resume completed work, and decide whether observable output
requirements have been met.

## Naming and layout

The skill follows the portable Agent Skills directory convention:

```text
run-zstar-workflows/
|-- SKILL.md
|-- agents/openai.yaml
|-- references/
`-- scripts/preflight.py
```

The folder and YAML `name` are both `run-zstar-workflows`: lowercase ASCII,
hyphen-separated, action-oriented, and shorter than 64 characters.

## Install

After installing ZStar, install the bundled skill into the default Codex skill
directory:

```bash
pip install -U zstar
zstar skill install
```

Open a new agent session so skill discovery is refreshed. To update an existing
copy after upgrading ZStar:

```bash
zstar skill install --force
```

Use a custom parent directory for another compatible harness:

```bash
zstar skill install --dest /path/to/skills
zstar skill path
```

The skill is also present in a source checkout at
`zstar/agent_skills/run-zstar-workflows/` and can be installed directly from
that repository path by an Agent Skill installer.

## Invoke

Explicit invocation uses the standard skill name:

```text
Use $run-zstar-workflows to inspect this BaTiO3 bulk BEC workspace, run a
preflight, generate a Slurm dry-run driver, and do not submit it.
```

The description also supports automatic discovery for requests that clearly
involve running or validating ZStar calculations.

## Machine-readable preflight

Agents should inspect a workspace before constructing or launching a workflow:

```bash
zstar skill preflight --root . --lane bec --dim bulk
zstar skill preflight --root . --lane raman --dim 2d
zstar skill preflight --root . --lane ir --dim molecule
zstar skill preflight --root . --lane database --dim 1d
```

The command writes JSON to standard output and does not modify the workspace.
It reports:

- ZStar and Python versions;
- selected scientific lane and dimensional convention;
- required input blockers and environment warnings;
- detected executables and retained artifacts;
- `.zstar/stages/*.json` status counts and failed-stage diagnostics; and
- a Boolean `ready` value based on required inputs, without pretending to
  certify physical convergence.

## Agent contract

The skill encodes the non-obvious invariants that should survive model changes:

| Contract | Agent behavior |
| --- | --- |
| Reference first | Complete and gate `0.no-move` before displaced BEC stages. |
| Dimensional split | Use Berry response along periodic axes and cube-integrated dipoles along open axes for 1D wires and 2D slabs. |
| 1D boundary | Run the implemented `z`-periodic ABACUS + PYATB BEC and Gamma-spectroscopy route, but do not claim finite-wavevector polar phonons without a genuine 1D Coulomb cutoff. |
| Resumability | Reuse `.zstar` state and repeat the same serial executor command. |
| Raman derivative | Use positive/negative normal-coordinate pairs. |
| Completion | Check named files, state records, units, and physical convention. |
| Authorization | Do not infer permission to submit or run expensive remote jobs. |

The result is a lightweight agent interface built from text instructions,
deterministic JSON inspection, stable CLI commands, and machine-verifiable
artifacts. It remains readable and usable without a particular language model.
