# Spectroscopy and MD lanes

## Infrared

For bulk or 2D systems:

```bash
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
  --dielectric BORN --dim 3 --outdir ir_spectrum
```

Use `--dim 2` for sheet response. A molecular IR calculation uses the same
positive/negative mode tree as Raman:

```bash
zstar ir --qpoints qpoints.yaml --dim 0 \
  --displacements raman --outdir ir_spectrum
```

## Raman

```bash
zstar raman prepare --stru STRU --qpoints qpoints.yaml \
  --modes "4-12" --outdir raman --copy INPUT-scf --copy KPT
zstar raman run --raman-dir raman --reference 0.no-move \
  --qpoints qpoints.yaml --dim 3 --dry-run
zstar raman status --raman-dir raman
```

After the dry run and environment check, repeat without `--dry-run`. Use
`--dim 2` for a sheet-susceptibility derivative or `--dim 0` for a molecule.
The molecular runner can write both channels with `--spectrum-outdir` and
`--ir-outdir`.

## VASP and CP2K calculators

For VASP, start from an `IBRION=5/6` vibrational `vasprun.xml` and converged
input directory:

```bash
zstar spectra prepare --calculator vasp --input-dir vasp_input \
  --modes-xml phonon/vasprun.xml --root vasp_spectra --dim 3
zstar spectra run --root vasp_spectra --command "mpirun -np 20 vasp_std" --dry-run
zstar spectra status --root vasp_spectra
```

For CP2K molecular native intensities:

```bash
zstar spectra prepare --calculator cp2k --input molecule.inp \
  --root cp2k_spectra --dim 0
zstar spectra run --root cp2k_spectra \
  --command "cp2k.ssmp -i input.inp -o output.log" --dry-run
```

Repeat without `--dry-run` only after inspecting the generated inputs. Finish
with `zstar spectra collect --root WORK`. Use `zstar spectra script` for one
shell, Slurm, or Torque driver. Native VASP/CP2K spectroscopy accepts `--dim 0`
or `--dim 3`, not `--dim 2`.

## MD + BEC

Use a fixed tensor set:

```bash
zstar md --dump dump.lammpstrj --fixed-bec Z-BORN-symm.out \
  --electronic-dielectric BORN --temperature 300 \
  --type-map "1:Hf,2:Zr,3:O" --outdir md_dielectric
```

Or use frame-resolved tensors:

```bash
zstar md --dump dump.lammpstrj --bec-dir bec_frames \
  --bec-pattern "frame_{step}.npy" --temperature 300 \
  --type-map "1:Hf,2:Zr,3:O" --outdir md_dielectric
```

Check equilibration, frame/step matching, atom order, periodic unwrapping, and
sampling convergence. Report ionic susceptibility and total dielectric response
separately. If no electronic tensor is supplied, identify the result as
`I + chi_ionic`, not a fully electronic-plus-ionic prediction.
