# zstar/cli.py
# -*- coding: utf-8 -*-

"""
Unified CLI for the ZStar toolkit.

Subcommands:
- gen, deal, born, polar, workflow
- ph, postph
- wyckoff, irrep, ir, raman, vasp
- calc, freq, md, potential
- symcheck, bornsym

This is adapted from the original pykappa-dev CLI, with imports rewritten to
use the installed zstar package modules and exposed as the `zstar` entry point.
"""

import argparse
import os
import sys

from . import __version__

VERSION_STR = f"ZStar {__version__}"


def zstar_cli(argv=None) -> None:
    """
    Entry point function for the `zstar` command.

    If `argv` is None, arguments are taken from `sys.argv[1:]` (normal CLI use).
    """
    parser = argparse.ArgumentParser(
        prog="zstar",
        description="ZStar: A Python toolkit for first-principles Born effective charge, "
                    "polarization, phonon, and dielectric analyses."
    )
    parser.add_argument('--version', action='store_true', help='Show version and exit')
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # ---------------- gen ----------------
    parser_gen = subparsers.add_parser('gen', help='Generate polarization data.')
    parser_gen.add_argument('-i', '--input', type=str, default=None,
                            help='Given your own INPUT file for ABACUS SCF')
    parser_gen.add_argument('--dim', type=int, help='Dim of your systems, 2 for 2D, default is 3.', default=3)
    parser_gen.add_argument('--method', type=str, help='Type of finite difference method, by forward or central, with pricesion of first and seconde order.', default='forward')
    parser_gen.add_argument('--xc', type=str,
                            help='dft_functional in abacus INPUT, default is pbe, you can change to pbesol',
                            default='pbe')
    parser_gen.add_argument('--vdw', type=str, help='vdw in INPUT', default=None)
    parser_gen.add_argument('--init', type=bool,
                            help='init_chg in INPUT, set False to use atomic',
                            default=True)
    parser_gen.add_argument('--kspacing', type=float,
                            help='kspacing in INPUT, default 0.1',
                            default=0.1)
    parser_gen.add_argument('--force', action='store_true',
                            help="Force overwrite existing directories if they exist.",
                            default=False)
    parser_gen.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_gen.add_argument('--symmprec', '--tol', type=float,
                            help='Symmetry precision of STRU, default is 1e-3',
                            default=1e-3)
    parser_gen.add_argument('--atom', type=str,
                            help='List of atoms separated by space',
                            default=None)
    parser_gen.add_argument('--move', type=str,
                            help='Directions (x, y, z) separated by no space',
                            default=None)
    parser_gen.add_argument('--input-mode',
                            choices=['abacus', 'pyatb', 'hamgnn', 'custom'],
                            default=None,
                            help='Input preparation mode; default follows calculator')
    parser_gen.add_argument('--input_sets',
                            default=None,
                            help='Extra input files or a directory. '
                                 'Space-separated list or an absolute directory path.')
    # 只保留其一：默认启用 reduce；若指定 --all 则关闭
    gen_scope = parser_gen.add_mutually_exclusive_group()
    gen_scope.add_argument(
        '--reduce', '--star',
        dest='reduce',
        action='store_true',
        default=True,
        help='Use only starred/reduced atoms (default).'
    )
    gen_scope.add_argument(
        '--all',
        dest='reduce',
        action='store_false',
        help='Use all atoms (disable reduction).'
    )
    gen_calc = parser_gen.add_mutually_exclusive_group()
    gen_calc.add_argument('--abacus', action='store_true',
                          help='Use ABACUS for NSCF Berry phase.')
    gen_calc.add_argument('--pyatb', action='store_true',
                          help='Use PyATB for NSCF Berry phase [Recommended].')

    # ---------------- deal ----------------
    parser_deal = subparsers.add_parser(
        'deal', help='Deal with polarization data to get BORN effective charge.'
    )
    parser_deal.add_argument('--dim', type=int,
                             help='Dim of your systems, 2 for 2D, default is 3.', default=3)
    parser_deal.add_argument('--method', type=str, help='Finite difference method, by forward or central, with pricesion of first and seconde order. To save calculation resource you can choose forward', default='forward')
    parser_deal.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_deal.add_argument('--symmprec', '--tol', type=float,
                             help='Symmetry precision of STRU, default is 1e-3',
                             default=1e-3)
    parser_deal.add_argument('--solo', action='store_true',
                             help='Running type of deal_polar: set solo for only polar')
    deal_calc = parser_deal.add_mutually_exclusive_group()
    deal_calc.add_argument('--abacus', action='store_true',
                           help='Use ABACUS for NSCF Berry phase.')
    deal_calc.add_argument('--pyatb', action='store_true',
                           help='Use PyATB for NSCF Berry phase [Recommended].')

    # ---------------- bornsym ----------------
    parser_borns = subparsers.add_parser(
        'bornsym',
        help='Generate full Born tensors from REDUCED via symmetry (no full reference required).'
    )
    parser_borns.add_argument('--stru', default='STRU',
                              help='Path to STRU used for symmetry (default: STRU)')
    parser_borns.add_argument('--reduced', default='Z-BORN-reduced.out',
                              help='Reduced Born file with starred atoms '
                                   '(default: Z-BORN-reduced.out)')
    parser_borns.add_argument('--symmprec', '--tol', type=float, default=1e-3,
                              help='Symmetry precision (default: 1e-3)')
    parser_borns.add_argument('--out', default='born_generation_from_symm.log',
                              help='Generation log output '
                                   '(default: born_generation_from_symm.log)')

    # ---------------- born ----------------
    parser_born = subparsers.add_parser(
        'born', help='Deal with polarization data to get BORN effective charge.'
    )
    parser_born.add_argument('--dim', type=int,
                             help='Dim of your systems, 2 for 2D.', default=3)
    parser_born.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_born.add_argument('--symmprec', '--tol', type=float,
                             help='Symmetry precision of STRU, default is 1e-3',
                             default=1e-3)
    parser_born.add_argument('--solo', action='store_true',
                             help='Running type of deal_polar: set solo for only polar')
    born_calc = parser_born.add_mutually_exclusive_group()
    born_calc.add_argument('--abacus', action='store_true',
                           help='Use ABACUS for NSCF Berry phase.')
    born_calc.add_argument('--pyatb', action='store_true',
                           help='Use PyATB for NSCF Berry phase [Recommended].')

    # ---------------- polar ----------------
    parser_polar = subparsers.add_parser('polar', help='Polarization data only.')
    parser_polar.add_argument('--dim', type=int,
                              help='Dim of your systems, 2 for 2D.', default=3)
    parser_polar.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_polar.add_argument('--symmprec', '--tol', type=float,
                              help='Symmetry precision of STRU, default is 1e-3',
                              default=1e-3)
    parser_polar.add_argument('--solo', action='store_true',
                              help='Running type of deal_polar: set solo for only polar')
    polar_calc = parser_polar.add_mutually_exclusive_group()
    polar_calc.add_argument('--abacus', action='store_true',
                            help='Use ABACUS for NSCF Berry phase.')
    polar_calc.add_argument('--pyatb', action='store_true',
                            help='Use PyATB for NSCF Berry phase [Recommended].')

    # ---------------- ph ----------------
    parser_ph = subparsers.add_parser('ph', help='Generate phonon data.')
    parser_ph.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_ph.add_argument('--symmprec', '--tol', type=float,
                           help='Symmetry precision of STRU, default is 1e-3',
                           default=1e-3)
    parser_ph.add_argument('--node', help='Node to sub', default='s1')
    parser_ph.add_argument('--dim', help='Dim of phonopy', default='1 1 1')

    # ---------------- postph ----------------
    parser_postph = subparsers.add_parser(
        'postph',
        help='Post-process phonon data. Auto detect phonopy_disp.yaml to get the DIM and SYMMPREC'
    )
    parser_postph.add_argument('--stru', help='Path to the STRU file', default='STRU')
    parser_postph.add_argument('--dim', help='Dim of phonopy', default=None)
    parser_postph.add_argument('--nac', action='store_true',
                               help='Whether to use NAC correction, default is False')
    parser_postph.add_argument('--symmprec', '--tol', type=float,
                               help='Symmetry precision of STRU, default is 1e-3',
                               default=1e-3)

    # ---------------- serial workflow ----------------
    parser_workflow = subparsers.add_parser(
        'workflow',
        help='Run or generate a serial, resumable Born-charge workflow.'
    )
    workflow_actions = parser_workflow.add_subparsers(
        dest='workflow_action', required=True
    )
    parser_workflow_run = workflow_actions.add_parser(
        'run', help='Run 0.no-move and all displacement stages serially.'
    )
    parser_workflow_run.add_argument('--root', default='.')
    parser_workflow_run.add_argument(
        '--abacus-command', default='mpirun -np 1 abacus'
    )
    parser_workflow_run.add_argument('--pyatb-input', default='pyatb_input')
    parser_workflow_run.add_argument(
        '--pyatb-command', default='mpirun -np 1 pyatb'
    )
    parser_workflow_run.add_argument('--pyatb-executable', default='pyatb')
    parser_workflow_run.add_argument('--mp-density', type=float, default=0.08)
    parser_workflow_run.add_argument(
        '--gap-mode', choices=['path', 'mp'], default='path',
        help='One-time 0.no-move PYATB band-path gate; use mp for stricter sampling.'
    )
    parser_workflow_run.add_argument(
        '--dimensionality', '--dim', type=int, choices=[2, 3], default=3
    )
    parser_workflow_run.add_argument('--min-gap', type=float, default=0.01)
    parser_workflow_run.add_argument(
        '--no-insulation-check', action='store_true',
        help='Disable the PYATB band-gap gate (not recommended).'
    )
    parser_workflow_run.add_argument(
        '--no-electronic-dielectric', action='store_true'
    )
    parser_workflow_run.add_argument(
        '--legacy-omega-max', type=float, default=30.0,
        help='Legacy PYATB optical cutoff in eV (default: 30; 0.1 eV spacing).'
    )
    parser_workflow_run.add_argument(
        '--legacy-domega', type=float, default=0.10
    )
    parser_workflow_run.add_argument('--omp-threads', type=int, default=1)
    parser_workflow_run.add_argument('--dry-run', action='store_true')
    parser_workflow_run.add_argument('--stop-after', type=int, default=None)

    parser_workflow_status = workflow_actions.add_parser(
        'status', help='Inspect stage progress and band-gap gates.'
    )
    parser_workflow_status.add_argument('--root', default='.')

    parser_workflow_script = workflow_actions.add_parser(
        'script', help='Generate one shell, Slurm, or Torque driver script.'
    )
    parser_workflow_script.add_argument(
        '--backend', choices=['shell', 'slurm', 'torque'], required=True
    )
    parser_workflow_script.add_argument('--root', default='.')
    parser_workflow_script.add_argument('--output', default=None)
    parser_workflow_script.add_argument('--job-name', default='zstar-born')
    parser_workflow_script.add_argument('--nodes', type=int, default=1)
    parser_workflow_script.add_argument('--tasks', type=int, default=1)
    parser_workflow_script.add_argument('--cpus-per-task', type=int, default=28)
    parser_workflow_script.add_argument('--walltime', default='24:00:00')
    parser_workflow_script.add_argument('--queue', default=None)
    parser_workflow_script.add_argument('--account', default=None)
    parser_workflow_script.add_argument('--env-script', default=None)
    parser_workflow_script.add_argument(
        '--abacus-command', default='mpirun -np 1 abacus'
    )
    parser_workflow_script.add_argument(
        '--pyatb-command', default='mpirun -np 1 pyatb'
    )
    parser_workflow_script.add_argument('--mp-density', type=float, default=0.08)
    parser_workflow_script.add_argument(
        '--gap-mode', choices=['path', 'mp'], default='path'
    )
    parser_workflow_script.add_argument(
        '--dimensionality', '--dim', type=int, choices=[2, 3], default=3
    )
    parser_workflow_script.add_argument('--min-gap', type=float, default=0.01)
    parser_workflow_script.add_argument(
        '--no-insulation-check', action='store_true'
    )
    parser_workflow_script.add_argument(
        '--legacy-omega-max', type=float, default=30.0
    )
    parser_workflow_script.add_argument(
        '--submit', action='store_true',
        help='Submit with sbatch/qsub, or run the shell script immediately.'
    )

    # ---------------- symcheck ----------------
    parser_symcheck = subparsers.add_parser(
        'symcheck',
        help='Verify Born tensors vs symmetry using a FULL reference (Z-BORN-all.out).'
    )
    parser_symcheck.add_argument('--stru', default='STRU',
                                 help='Path to STRU used for symmetry (default: STRU)')
    parser_symcheck.add_argument('--reduced', default='Z-BORN-reduced.out',
                                 help='Reduced Born file with starred atoms '
                                      '(default: Z-BORN-reduced.out)')
    parser_symcheck.add_argument('--allfile', default='Z-BORN-all.out', required=False,
                                 help='Full reference Born tensor file '
                                      '(e.g., Z-BORN-all.out)')
    parser_symcheck.add_argument('--symmprec', '--tol', type=float, default=1e-3,
                                 help='Symmetry precision (default: 1e-3)')
    parser_symcheck.add_argument('--out', default='born_symmetry_report.txt',
                                 help='Text report output '
                                      '(default: born_symmetry_report.txt)')
    parser_symcheck.add_argument('--json', dest='json_path', default='born_symmetry_report.json',
                                 help='JSON report output '
                                      '(default: born_symmetry_report.json)')
    parser_symcheck.add_argument('--csv', dest='csv_path', default=None,
                                 help='Optional CSV dump path')

    # ---------------- wyckoff / irrep / vasp / calc ----------------
    parser_wyckoff = subparsers.add_parser('wyckoff', help='Get Wyckoff positions.')
    parser_wyckoff.add_argument('--stru', help='Path to the STRU file', default='STRU')

    p_ir = subparsers.add_parser('irrep', help='Classify Gamma irreps from irreps.yaml')
    p_ir.add_argument('-f', '--file', default='irreps.yaml', help='Path to irreps.yaml')
    p_ir.add_argument('--mode', default='db', choices=['db', 'default', 'smodes'],
                      help='Classification mode: db/default (no external tools) or '
                           'smodes (legacy)')
    p_ir.add_argument('--acoustic-thz', type=float, default=0.05,
                      help='Acoustic threshold (THz)')
    p_ir.add_argument('--stru', default=None,
                      help='(compat) STRU path; smodes will look in CWD anyway')

    parser_ir = subparsers.add_parser(
        'ir',
        help='Calculate mode effective charges and infrared spectra.'
    )
    parser_ir.add_argument('--qpoints', default='qpoints.yaml')
    parser_ir.add_argument(
        '--born', default='Z-BORN-symm.out',
        help='Z-BORN-symm.out or a Phonopy-style BORN file.'
    )
    parser_ir.add_argument(
        '--dielectric', default=None,
        help='Optional BORN or PYATB optical output used for epsilon infinity.'
    )
    parser_ir.add_argument('--dim', type=int, choices=[2, 3], default=3)
    parser_ir.add_argument(
        '--thickness', type=float, default=None,
        help='2D effective thickness in Angstrom; omit for sheet polarizability.'
    )
    parser_ir.add_argument(
        '--modes', default=None,
        help='One-based mode list, e.g. 4,5,8-10. Default: all optical modes.'
    )
    parser_ir.add_argument('--acoustic-cutoff', type=float, default=5.0)
    parser_ir.add_argument('--broadening', type=float, default=10.0)
    parser_ir.add_argument('--max-frequency', type=float, default=None)
    parser_ir.add_argument('--points', type=int, default=2001)
    parser_ir.add_argument('--outdir', default='ir_spectrum')
    parser_ir.add_argument('--no-plot', action='store_true')

    parser_raman = subparsers.add_parser(
        'raman',
        help='Prepare finite differences or calculate Raman spectra.'
    )
    raman_actions = parser_raman.add_subparsers(
        dest='raman_action', required=True
    )
    parser_raman_prepare = raman_actions.add_parser(
        'prepare', help='Generate +/- Gamma-mode displacement structures.'
    )
    parser_raman_prepare.add_argument('--stru', default='STRU')
    parser_raman_prepare.add_argument('--qpoints', default='qpoints.yaml')
    parser_raman_prepare.add_argument('--outdir', default='raman')
    parser_raman_prepare.add_argument(
        '--amplitude', type=float, default=0.02,
        help='Normal-coordinate step in Angstrom sqrt(amu).'
    )
    parser_raman_prepare.add_argument('--modes', default=None)
    parser_raman_prepare.add_argument('--acoustic-cutoff', type=float, default=5.0)
    parser_raman_prepare.add_argument(
        '--copy', action='append', default=None,
        help='Additional calculation input copied into every +/- directory.'
    )

    parser_raman_run = raman_actions.add_parser(
        'run',
        help='Run all +/- structures serially, collect tensors, and make a spectrum.'
    )
    parser_raman_run.add_argument('--raman-dir', default='raman')
    parser_raman_run.add_argument('--reference', default='0.no-move')
    parser_raman_run.add_argument('--qpoints', default='qpoints.yaml')
    parser_raman_run.add_argument('--dim', type=int, choices=[2, 3], default=3)
    parser_raman_run.add_argument(
        '--abacus-command', default='mpirun -np 1 abacus'
    )
    parser_raman_run.add_argument('--pyatb-input', default='pyatb_input')
    parser_raman_run.add_argument(
        '--pyatb-command', default='mpirun -np 1 pyatb'
    )
    parser_raman_run.add_argument('--pyatb-executable', default='pyatb')
    parser_raman_run.add_argument('--mp-density', type=float, default=0.08)
    parser_raman_run.add_argument(
        '--gap-mode', choices=['path', 'mp'], default='path',
        help='One-time reference PYATB band-path gate; use mp for stricter sampling.'
    )
    parser_raman_run.add_argument('--min-gap', type=float, default=0.01)
    parser_raman_run.add_argument(
        '--no-insulation-check', action='store_true',
        help='Disable the PYATB band-gap gate (not recommended).'
    )
    parser_raman_run.add_argument('--legacy-omega-max', type=float, default=30.0)
    parser_raman_run.add_argument('--legacy-domega', type=float, default=0.10)
    parser_raman_run.add_argument('--omp-threads', type=int, default=1)
    parser_raman_run.add_argument('--dry-run', action='store_true')
    parser_raman_run.add_argument('--stop-after', type=int, default=None)
    parser_raman_run.add_argument('--temperature', type=float, default=300.0)
    parser_raman_run.add_argument('--laser-nm', type=float, default=532.0)
    parser_raman_run.add_argument('--broadening', type=float, default=8.0)
    parser_raman_run.add_argument('--max-frequency', type=float, default=None)
    parser_raman_run.add_argument('--points', type=int, default=2001)
    parser_raman_run.add_argument(
        '--spectrum-outdir', default='raman_spectrum'
    )
    parser_raman_run.add_argument('--no-spectrum', action='store_true')
    parser_raman_run.add_argument('--no-plot', action='store_true')

    parser_raman_status = raman_actions.add_parser(
        'status', help='Inspect Raman +/- stage progress.'
    )
    parser_raman_status.add_argument('--raman-dir', default='raman')

    parser_raman_collect = raman_actions.add_parser(
        'collect', help='Central-difference PYATB dielectric tensors.'
    )
    parser_raman_collect.add_argument('--raman-dir', default='raman')
    parser_raman_collect.add_argument('--qpoints', default='qpoints.yaml')
    parser_raman_collect.add_argument(
        '--dim', type=int, choices=[2, 3], default=3
    )

    parser_raman_spectrum = raman_actions.add_parser(
        'spectrum', help='Calculate Placzek Raman activities and spectrum.'
    )
    parser_raman_spectrum.add_argument('--qpoints', default='qpoints.yaml')
    raman_tensor_source = parser_raman_spectrum.add_mutually_exclusive_group(
        required=True
    )
    raman_tensor_source.add_argument('--tensors', default=None)
    raman_tensor_source.add_argument('--raman-dir', default=None)
    parser_raman_spectrum.add_argument(
        '--dim', type=int, choices=[2, 3], default=3
    )
    parser_raman_spectrum.add_argument('--temperature', type=float, default=300.0)
    parser_raman_spectrum.add_argument('--laser-nm', type=float, default=532.0)
    parser_raman_spectrum.add_argument('--broadening', type=float, default=8.0)
    parser_raman_spectrum.add_argument('--max-frequency', type=float, default=None)
    parser_raman_spectrum.add_argument('--points', type=int, default=2001)
    parser_raman_spectrum.add_argument('--outdir', default='raman_spectrum')
    parser_raman_spectrum.add_argument('--no-plot', action='store_true')

    parser_vasp = subparsers.add_parser(
        'vasp', help='Convert ABACUS structure format STRU to VASP format POSCAR.'
    )
    parser_vasp.add_argument('--stru', help='Path to the STRU file', default='STRU')

    parser_calc = subparsers.add_parser('calc', help='Calculate static dielectric tensor.')
    parser_calc.add_argument('--tolerance', type=float,
                             help="Set the tolerance of ZERO in dielectric tensor matrix.",
                             default=1e-3)
    parser_calc.add_argument('--ir-tolerance', type=float,
                             help="Set the tolerance of infrared modes of all modes.",
                             default=5e-2)
    parser_calc.add_argument('--ir-choose',
                             help="Set choose only infrared modes or all modes.",
                             choices=['ir', 'all'], default='ir')
    parser_calc.add_argument('--plot', action='store_true',
                             help="Whether to plot the dielectric constant VS the frequency.",
                             default=False)
    # new flags for updated pipeline
    parser_calc.add_argument('--mode', default='db',
                             choices=['db', 'default', 'smodes'],
                             help="Mode for classification (default: db).")
    parser_calc.add_argument('--stru', dest='stru_file', default='STRU',
                             help='Structure file (STRU or POSCAR/vasp).')
    parser_calc.add_argument('--irreps', dest='irreps_file', default='irreps.yaml',
                             help='Path to irreps.yaml.')
    parser_calc.add_argument('--qpoints', default='qpoints.yaml')
    parser_calc.add_argument('--born', default='BORN')
    parser_calc.add_argument('--dielectric', default=None)
    parser_calc.add_argument('--dim', type=int, choices=[2, 3], default=3)
    parser_calc.add_argument('--thickness', type=float, default=None)
    parser_calc.add_argument('--acoustic-cutoff', type=float, default=5.0)
    parser_calc.add_argument('--broadening', type=float, default=10.0)
    parser_calc.add_argument('--max-frequency', type=float, default=None)
    parser_calc.add_argument('--points', type=int, default=2001)
    parser_calc.add_argument('--outdir', default='dielectric_response')

    parser_freq = subparsers.add_parser('freq', help='Calculate the frequency-dependent dielectric functions.')
    parser_freq.add_argument('--tolerance', type=float,
                             help="Set the tolerance of ZERO in dielectric tensor matrix.",
                             default=1e-3)
    parser_freq.add_argument('--ir-tolerance', type=float,
                             help="Set the tolerance of infrared modes of all modes.",
                             default=5e-2)
    parser_freq.add_argument('--ir-choose',
                             help="Set choose only infrared modes or all modes.",
                             choices=['ir', 'all'], default='ir')
    parser_freq.add_argument('--plot', action='store_true',
                             help="Whether to plot the dielectric constant VS the frequency.",
                             default=True)
    # new flags for updated pipeline
    parser_freq.add_argument('--mode', default='db',
                             choices=['db', 'default', 'smodes'],
                             help="Mode for classification (default: db).")
    parser_freq.add_argument('--stru', dest='stru_file', default='STRU',
                             help='Structure file (STRU or POSCAR/vasp).')
    parser_freq.add_argument('--irreps', dest='irreps_file', default='irreps.yaml',
                             help='Path to irreps.yaml.')
    parser_freq.add_argument('--qpoints', default='qpoints.yaml')
    parser_freq.add_argument('--born', default='BORN')
    parser_freq.add_argument('--dielectric', default=None)
    parser_freq.add_argument('--dim', type=int, choices=[2, 3], default=3)
    parser_freq.add_argument('--thickness', type=float, default=None)
    parser_freq.add_argument('--acoustic-cutoff', type=float, default=5.0)
    parser_freq.add_argument('--broadening', type=float, default=10.0)
    parser_freq.add_argument('--max-frequency', type=float, default=None)
    parser_freq.add_argument('--points', type=int, default=2001)
    parser_freq.add_argument('--outdir', default='dielectric_response')

    parser_md = subparsers.add_parser(
        'md',
        help='Calculate total dielectric tensor from MD trajectory and supplied BEC tensors.'
    )
    md_source = parser_md.add_mutually_exclusive_group(required=True)
    md_source.add_argument('--dump', '--dump-file', dest='dump_file',
                           help='LAMMPS dump trajectory.')
    md_source.add_argument('--structure-dir',
                           help='Directory containing structure frames such as frame_*.vasp.')
    md_bec = parser_md.add_mutually_exclusive_group(required=True)
    md_bec.add_argument('--bec-dir',
                        help='Directory containing per-frame BEC tensors.')
    md_bec.add_argument('--fixed-bec',
                        help='Fixed BEC tensor file used for all frames.')
    parser_md.add_argument('--bec-pattern', default='frame_{step}.npy',
                           help='Per-frame BEC filename pattern.')
    parser_md.add_argument('--structure-glob', default='frame_*.vasp',
                           help='Glob pattern for structure frames.')
    parser_md.add_argument('--temperature', type=float, required=True,
                           help='MD temperature in K.')
    parser_md.add_argument('--type-map', default='',
                           help="LAMMPS type map, e.g. '1:Hf,2:Zr,3:O'.")
    parser_md.add_argument('--start-step', type=int, default=None,
                           help='First MD step to include.')
    parser_md.add_argument('--end-step', type=int, default=None,
                           help='Last MD step to include.')
    parser_md.add_argument('--stride-step', type=int, default=1,
                           help='Step-space stride after range filtering.')
    parser_md.add_argument('--second-half', action='store_true',
                           help='Use only the second half of selected frames.')
    parser_md.add_argument('--reference', choices=['mean', 'first'], default='mean',
                           help='Reference structure for displacements.')
    parser_md.add_argument('--remove-global-translation', action='store_true',
                           help='Remove per-frame centroid motion before computing displacements.')
    parser_md.add_argument('--no-minimum-image', action='store_true',
                           help='Disable minimum-image displacement reconstruction.')
    parser_md.add_argument('--volume', dest='volume_A3', type=float, default=None,
                           help='Override average volume in Angstrom^3.')
    parser_md.add_argument('--max-step-gap', type=int, default=None,
                           help='Maximum allowed step gap when matching nearest BEC frame.')
    parser_md.add_argument('--raw-moment-average', action='store_true',
                           help='Use <M M^T> instead of fluctuation covariance.')
    parser_md.add_argument('--unbiased', action='store_true',
                           help='Use N-1 covariance normalization.')
    parser_md.add_argument(
        '--electronic-dielectric', '--epsilon-infinity',
        dest='electronic_dielectric',
        default=None,
        help='BORN file or PYATB output containing epsilon infinity.'
    )
    parser_md.add_argument('--outdir', default='md_dielectric',
                           help='Output directory.')

    parser_potential = subparsers.add_parser(
        'potential',
        aliases=['pot'],
        help='Analyze ABACUS electrostatic-potential cube files.'
    )
    parser_potential.add_argument('--cube', default=None,
                                  help='Path to ElecStaticPot.cube. Auto-detects by default.')
    parser_potential.add_argument('--outdir', default=None,
                                  help='Output directory. Defaults to the cube directory.')
    parser_potential.add_argument('--prefix', default='ElecStaticPot',
                                  help='Output filename prefix.')
    parser_potential.add_argument('--axis', action='append',
                                  choices=['x', 'y', 'z'],
                                  help='Axis for plane-averaged 1D profile.')
    parser_potential.add_argument('--axes', default=None,
                                  help='Compact axis list, e.g. xyz or z.')
    parser_potential.add_argument('--plane', action='append',
                                  choices=['xy', 'xz', 'yz'],
                                  help='Plane map to write.')
    parser_potential.add_argument('--plane-index', type=int, default=None,
                                  help='Grid index along the plane normal.')
    parser_potential.add_argument('--plane-fraction', type=float, default=None,
                                  help='Fractional normal position from 0 to 1.')
    parser_potential.add_argument('--plane-average', action='store_true',
                                  help='Average the potential over the plane normal.')
    parser_potential.add_argument('--plane-coords',
                                  choices=['cartesian', 'rectified'],
                                  default='cartesian',
                                  help='Coordinate system for plane maps.')
    parser_potential.add_argument('--tile', nargs=2, type=int,
                                  metavar=('NA', 'NB'), default=(1, 1),
                                  help='Tile plane-map plots by NA x NB cells, e.g. --tile 5 5.')
    parser_potential.add_argument('--no-cell-frame', action='store_true',
                                  help='Do not draw the dashed central-unit-cell frame on tiled maps.')
    parser_potential.add_argument('--direction', action='append', default=None,
                                  help='Lattice direction for a perpendicular-plane-averaged profile, e.g. a+b or a-b.')
    parser_potential.add_argument('--direction-bins', type=int, default=None,
                                  help='Number of bins used for each directional profile.')
    parser_potential.add_argument('--direction-tile-radius', type=int, default=1,
                                  help='Neighbor-cell radius used by the legacy binning method.')
    parser_potential.add_argument('--direction-method', action='append', default=None,
                                  choices=['bin', 'nearest', 'linear', 'cubic', 'all'],
                                  help='Directional averaging method. Repeat to compare methods; default is linear.')
    parser_potential.add_argument('--direction-samples', nargs=2, type=int,
                                  metavar=('NU', 'NV'), default=(64, 64),
                                  help='Perpendicular-plane sample grid for interpolated directional profiles.')
    parser_potential.add_argument('--direction-smooth', type=float, default=0.0,
                                  help='Optional periodic Gaussian smoothing sigma in Angstrom for directional profiles.')
    parser_potential.add_argument('--value-unit',
                                  choices=['ry', 'ev', 'hartree'],
                                  default='ry',
                                  help='Potential unit stored in the cube data.')
    parser_potential.add_argument('--length-unit',
                                  choices=['bohr', 'angstrom'],
                                  default='bohr',
                                  help='Length unit used by cube header coordinates.')
    parser_potential.add_argument('--vacuum-level', action='store_true',
                                  help='Estimate z-direction vacuum level from the largest atom-free gap.')
    parser_potential.add_argument('--vacuum-sides', action='store_true',
                                  help='Estimate lower/upper z-vacuum levels and their potential step.')
    parser_potential.add_argument('--vacuum-exclude', type=float, default=6.0,
                                  help='Distance in Angstrom excluded from both sides of the vacuum gap.')
    parser_potential.add_argument('--center-slab', nargs='?', const='z',
                                  choices=['x', 'y', 'z'], default=None,
                                  help='Periodically shift the slab so its atomic center lies at the cell center.')
    parser_potential.add_argument('--polar-arrow',
                                  choices=['none', 'auto', '+x', '-x', '+y', '-y', '+z', '-z'],
                                  default='none',
                                  help='Draw a polarization-direction arrow on compatible profile plots.')
    parser_potential.add_argument('--no-plot', action='store_true',
                                  help='Write data files only.')
    parser_potential.add_argument('--dpi', type=int, default=300,
                                  help='Plot resolution.')
    parser_potential.add_argument('--cmap', default='viridis',
                                  help='Matplotlib colormap for plane maps.')

    args = parser.parse_args(argv)

    if args.version:
        print(VERSION_STR)
        return

    def _build_irrep_argv(a):
        """把 irrep 子命令解析到的参数转换成 read_irrep.main(argv) 需要的 argv 列表。"""
        argv_ir = []
        if getattr(a, 'file', None):
            argv_ir += ['--file', a.file]
        if getattr(a, 'mode', None):
            argv_ir += ['--mode', a.mode]
        if getattr(a, 'acoustic_thz', None) is not None:
            argv_ir += ['--acoustic-thz', str(a.acoustic_thz)]
        if getattr(a, 'stru', None):
            argv_ir += ['--stru', a.stru]
        return argv_ir or None  # 为空时传 None 让其自行处理

    def _parse_modes(text):
        if text is None or not str(text).strip():
            return None
        values = []
        for token in str(text).replace(' ', '').split(','):
            if not token:
                continue
            if '-' in token:
                start, stop = (int(value) for value in token.split('-', 1))
                values.extend(range(start, stop + 1))
            else:
                values.append(int(token))
        return list(dict.fromkeys(values))

    # ---------------- dispatch ----------------
    if args.command == 'workflow':
        from .workflow import (
            format_status_table,
            generate_backend_script,
            run_serial_workflow,
            submit_backend_script,
            workflow_status,
        )

        if args.workflow_action == 'run':
            states = run_serial_workflow(
                root=args.root,
                abacus_command=args.abacus_command,
                pyatb_input=args.pyatb_input,
                pyatb_command=args.pyatb_command,
                pyatb_executable=args.pyatb_executable,
                mp_density=args.mp_density,
                electronic_dielectric=not args.no_electronic_dielectric,
                check_insulating=not args.no_insulation_check,
                gap_mode=args.gap_mode,
                dimensionality=args.dimensionality,
                min_gap_eV=args.min_gap,
                legacy_omega_max=args.legacy_omega_max,
                legacy_domega=args.legacy_domega,
                omp_threads=args.omp_threads,
                dry_run=args.dry_run,
                stop_after=args.stop_after,
            )
            print(format_status_table(states))
        elif args.workflow_action == 'status':
            print(format_status_table(workflow_status(args.root)))
        elif args.workflow_action == 'script':
            script = generate_backend_script(
                root=args.root,
                backend=args.backend,
                output=args.output,
                job_name=args.job_name,
                nodes=args.nodes,
                tasks=args.tasks,
                cpus_per_task=args.cpus_per_task,
                walltime=args.walltime,
                queue=args.queue,
                account=args.account,
                env_script=args.env_script,
                abacus_command=args.abacus_command,
                pyatb_command=args.pyatb_command,
                mp_density=args.mp_density,
                check_insulating=not args.no_insulation_check,
                gap_mode=args.gap_mode,
                dimensionality=args.dimensionality,
                min_gap_eV=args.min_gap,
                legacy_omega_max=args.legacy_omega_max,
            )
            print(f"[OUT] {script}")
            if args.submit:
                print(submit_backend_script(script, args.backend))

    elif args.command == 'ir':
        from .spectra import (
            calculate_ir_spectrum,
            load_gamma_modes,
            read_born_data,
            write_ir_outputs,
        )

        modes = load_gamma_modes(args.qpoints)
        born = read_born_data(
            args.born,
            natoms=len(modes.masses_amu),
            dielectric_path=args.dielectric,
        )
        result = calculate_ir_spectrum(
            modes,
            born,
            dimensionality=args.dim,
            mode_numbers=_parse_modes(args.modes),
            acoustic_cutoff_cm1=args.acoustic_cutoff,
            broadening_cm1=args.broadening,
            max_frequency_cm1=args.max_frequency,
            points=args.points,
            thickness_angstrom=args.thickness,
        )
        summary = write_ir_outputs(
            args.outdir, result, plot=not args.no_plot
        )
        print(
            f"Calculated {summary['modes']} IR modes; "
            f"response: {summary['response_kind']}"
        )
        print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command == 'raman':
        from .spectra import (
            calculate_raman_spectrum,
            collect_raman_tensors,
            load_gamma_modes,
            load_raman_tensors,
            prepare_raman_displacements,
            write_raman_outputs,
        )

        if args.raman_action == 'status':
            from .workflow import format_status_table, raman_workflow_status

            print(format_status_table(raman_workflow_status(args.raman_dir)))
            return

        modes = load_gamma_modes(args.qpoints)
        if args.raman_action == 'prepare':
            manifest = prepare_raman_displacements(
                args.stru,
                modes,
                args.outdir,
                amplitude=args.amplitude,
                mode_numbers=_parse_modes(args.modes),
                acoustic_cutoff_cm1=args.acoustic_cutoff,
                copy_files=args.copy,
            )
            print(f"Generated {len(manifest['modes'])} Raman mode pairs.")
            print(f"[OUT] {os.path.abspath(args.outdir)}")
        elif args.raman_action == 'run':
            from .workflow import format_status_table, run_raman_workflow

            states = run_raman_workflow(
                args.raman_dir,
                reference_dir=args.reference,
                abacus_command=args.abacus_command,
                pyatb_input=args.pyatb_input,
                pyatb_command=args.pyatb_command,
                pyatb_executable=args.pyatb_executable,
                mp_density=args.mp_density,
                check_insulating=not args.no_insulation_check,
                gap_mode=args.gap_mode,
                dimensionality=args.dim,
                min_gap_eV=args.min_gap,
                legacy_omega_max=args.legacy_omega_max,
                legacy_domega=args.legacy_domega,
                omp_threads=args.omp_threads,
                dry_run=args.dry_run,
                stop_after=args.stop_after,
            )
            print(format_status_table(states))
            if args.dry_run or args.stop_after is not None:
                return
            mode_numbers, tensors, tensor_kind = collect_raman_tensors(
                args.raman_dir,
                dimensionality=args.dim,
                cell_height_angstrom=modes.cell_height_angstrom,
            )
            print(f"Collected {len(mode_numbers)} {tensor_kind} tensors.")
            if not args.no_spectrum:
                result = calculate_raman_spectrum(
                    modes,
                    mode_numbers,
                    tensors,
                    tensor_kind=tensor_kind,
                    temperature_K=args.temperature,
                    laser_nm=args.laser_nm,
                    broadening_cm1=args.broadening,
                    max_frequency_cm1=args.max_frequency,
                    points=args.points,
                )
                summary = write_raman_outputs(
                    args.spectrum_outdir,
                    result,
                    plot=not args.no_plot,
                )
                print(f"Calculated {summary['modes']} Raman modes.")
                print(f"[OUT] {os.path.abspath(args.spectrum_outdir)}")
        elif args.raman_action == 'collect':
            mode_numbers, tensors, tensor_kind = collect_raman_tensors(
                args.raman_dir,
                dimensionality=args.dim,
                cell_height_angstrom=modes.cell_height_angstrom,
            )
            print(
                f"Collected {len(mode_numbers)} {tensor_kind} tensors "
                f"with shape {tensors.shape}."
            )
            print(f"[OUT] {os.path.abspath(args.raman_dir)}")
        elif args.raman_action == 'spectrum':
            if args.raman_dir:
                mode_numbers, tensors, tensor_kind = collect_raman_tensors(
                    args.raman_dir,
                    dimensionality=args.dim,
                    cell_height_angstrom=modes.cell_height_angstrom,
                )
            else:
                mode_numbers, tensors, tensor_kind = load_raman_tensors(
                    args.tensors
                )
            result = calculate_raman_spectrum(
                modes,
                mode_numbers,
                tensors,
                tensor_kind=tensor_kind,
                temperature_K=args.temperature,
                laser_nm=args.laser_nm,
                broadening_cm1=args.broadening,
                max_frequency_cm1=args.max_frequency,
                points=args.points,
            )
            summary = write_raman_outputs(
                args.outdir, result, plot=not args.no_plot
            )
            print(f"Calculated {summary['modes']} Raman modes.")
            print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command == 'irrep':
        from .read_irrep import main as run_read_irrep_cli

        run_read_irrep_cli(_build_irrep_argv(args))

    elif args.command == 'calc':
        from .spectra import (
            calculate_ir_spectrum,
            load_gamma_modes,
            read_born_data,
            write_ir_outputs,
        )

        modes = load_gamma_modes(args.qpoints)
        born = read_born_data(
            args.born,
            natoms=len(modes.masses_amu),
            dielectric_path=args.dielectric,
        )
        result = calculate_ir_spectrum(
            modes,
            born,
            dimensionality=args.dim,
            acoustic_cutoff_cm1=args.acoustic_cutoff,
            broadening_cm1=args.broadening,
            max_frequency_cm1=args.max_frequency,
            points=args.points,
            thickness_angstrom=args.thickness,
        )
        write_ir_outputs(args.outdir, result, plot=args.plot)
        print(f"Static {result.response_kind}:")
        for row in result.response_real[0]:
            print("  " + " ".join(f"{value:14.8e}" for value in row))
        print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command == 'freq':
        from .spectra import (
            calculate_ir_spectrum,
            load_gamma_modes,
            read_born_data,
            write_ir_outputs,
        )
        modes = load_gamma_modes(args.qpoints)
        born = read_born_data(
            args.born,
            natoms=len(modes.masses_amu),
            dielectric_path=args.dielectric,
        )
        result = calculate_ir_spectrum(
            modes,
            born,
            dimensionality=args.dim,
            acoustic_cutoff_cm1=args.acoustic_cutoff,
            broadening_cm1=args.broadening,
            max_frequency_cm1=args.max_frequency,
            points=args.points,
            thickness_angstrom=args.thickness,
        )
        write_ir_outputs(args.outdir, result, plot=True)
        print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command == 'md':
        from .md_dielectric import compute_md_dielectric, parse_type_map

        result = compute_md_dielectric(
            dump_file=args.dump_file,
            structure_dir=args.structure_dir,
            structure_glob=args.structure_glob,
            bec_dir=args.bec_dir,
            bec_pattern=args.bec_pattern,
            fixed_bec=args.fixed_bec,
            temperature=args.temperature,
            type_map=parse_type_map(args.type_map),
            start_step=args.start_step,
            end_step=args.end_step,
            stride_step=args.stride_step,
            second_half=args.second_half,
            reference=args.reference,
            remove_global_translation=args.remove_global_translation,
            minimum_image=not args.no_minimum_image,
            volume_A3=args.volume_A3,
            max_step_gap=args.max_step_gap,
            raw_moment_average=args.raw_moment_average,
            unbiased=args.unbiased,
            electronic_dielectric=args.electronic_dielectric,
            outdir=args.outdir,
        )
        print("Total dielectric tensor from MD dipole fluctuations:")
        for row in result.epsilon:
            print("  " + " ".join(f"{value:14.8e}" for value in row))
        print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command in ('potential', 'pot'):
        from .potential import analyze_potential, normalize_axes

        summary = analyze_potential(
            cube=args.cube,
            outdir=args.outdir,
            prefix=args.prefix,
            axes=normalize_axes(args.axis, args.axes),
            planes=args.plane,
            plane_index=args.plane_index,
            plane_fraction=args.plane_fraction,
            plane_average=args.plane_average,
            plane_coord_mode=args.plane_coords,
            tile=args.tile,
            highlight_cell=not args.no_cell_frame,
            directions=args.direction,
            direction_bins=args.direction_bins,
            direction_tile_radius=args.direction_tile_radius,
            direction_methods=args.direction_method,
            direction_samples=args.direction_samples,
            direction_smooth=args.direction_smooth,
            value_unit=args.value_unit,
            length_unit=args.length_unit,
            vacuum_level=args.vacuum_level,
            vacuum_sides=args.vacuum_sides,
            vacuum_exclude=args.vacuum_exclude,
            center_slab_axis=args.center_slab,
            polar_arrow=args.polar_arrow,
            plot=not args.no_plot,
            dpi=args.dpi,
            cmap=args.cmap,
        )
        print(f"Processed cube: {summary['cube']}")
        print(f"[OUT] {summary['outdir']}")
        for axis, info in summary["axis_profiles"].items():
            print(f"[AXIS {axis.upper()}] {info['dat']}")
        for plane, info in summary["plane_maps"].items():
            print(f"[PLANE {plane.upper()}] {info['dat']}")
        for direction, info in summary["direction_profiles"].items():
            target = info.get('dat') or info.get('png')
            print(f"[DIRECTION {direction}] {target}")

    elif args.command == 'gen':
        from .gen_polar import gen_polar as run_gen

        # Full x/y/z displacements are required for the hybrid 2D BEC tensor.
        if args.move is None or str(args.move).strip() == "":
            args.move = "x y z"
            print(
                f"[INFO] dim={args.dim} and --move not specified; "
                "default to full x/y/z displacements."
            )

        move_input = [c for c in str(args.move) if c in ('x', 'y', 'z')]
        print("处理后的 --move 参数:", move_input)

        nscf_calculator = 'abacus' if getattr(args, 'abacus', False) else (
                          'pyatb' if getattr(args, 'pyatb', False) else 'pyatb')
        input_mode = args.input_mode or nscf_calculator
        input_sets = args.input_sets

        method = getattr(args, "method", "forward")
        if method in ["center", "central"]:
            method_fd = "central"
        else:
            method_fd = "forward"

        run_gen(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            force_delete=args.force,
            atom_input=args.atom,
            move_input=move_input,
            scf_input=args.input,
            xc=args.xc,
            dimension=args.dim,
            vdw=args.vdw,
            init_chg_bool=args.init,
            k_grid=args.kspacing,
            nscf_calculator=nscf_calculator,
            input_mode=input_mode,
            input_sets=input_sets,
            extract_starred_atoms_only=args.reduce,
            method=method_fd
        )

    elif args.command == 'ph':
        from .phonon_gen import run_phonopy_and_process_files as run_ph

        run_ph(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            dim=args.dim,
            abacus_sub="abacus_x.sh",
            vasp_sub="vasp_scf.sh",
            node=args.node
        )

    elif args.command == 'postph':
        from .phonon_post import run_eigen_irrep as run_postph

        run_postph(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            nac=args.nac,
            dim=args.dim
        )

    elif args.command in ('rpolar', 'deal', 'born', 'polar'):
        from .deal_polar import main as run_deal_polar

        calc_flag = 'abacus' if getattr(args, 'abacus', False) else 'pyatb'
        running_type = 'solo' if getattr(args, 'solo', False) else None

        method = getattr(args, "method", "forward")
        if method in ["center", "central"]:
            method_fd = "central"
        else:
            method_fd = "forward"

        kwargs = dict(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            dimension=args.dim,
            method=method_fd,
            running_type=running_type
        )
        if calc_flag:
            kwargs['nscf_calculator'] = calc_flag

        run_deal_polar(**kwargs)

    elif args.command == 'vasp':
        from .get_wyckoff import stru2vasp as run_stru2vasp

        run_stru2vasp(f_stru=args.stru)

    elif args.command == 'wyckoff':
        from .get_wyckoff import get_wyckoff_position as run_get_wyckoff

        run_get_wyckoff(fstru=args.stru)

    elif args.command == 'symcheck':
        from .verify_born_symmetry import run_symcheck

        if not os.path.isfile(args.allfile):
            print(
                f"[ERROR] --allfile not found: {args.allfile}. "
                f"Please provide a full-atom Born file (e.g., Z-BORN-all.out).",
                file=sys.stderr
            )
            sys.exit(2)
        run_symcheck(
            stru=args.stru,
            reduced=args.reduced,
            all=args.allfile,
            symprec=args.symmprec,
            out=args.out,
            json_path=args.json_path,
            csv_path=args.csv_path
        )

    elif args.command == 'bornsym':
        from .verify_born_symmetry import run_symcheck

        # reduced-only generation (no full reference)
        run_symcheck(
            stru=args.stru,
            reduced=args.reduced,
            all=None,  # generation mode
            symprec=args.symmprec,
            out=args.out,
            csv_path=None
        )

    else:
        parser.print_help()


# 可选：方便你在源码树里直接 python -m zstar.cli 调试
if __name__ == "__main__":
    zstar_cli()
