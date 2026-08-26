# zstar/cli.py
# -*- coding: utf-8 -*-

"""
Unified CLI for the ZStar toolkit.

Subcommands:
- gen, deal, born, polar, polar2d, workflow
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


def _add_qe_action_parsers(action_subparsers) -> None:
    """Register the shared Quantum ESPRESSO response actions."""
    parser_prepare = action_subparsers.add_parser('prepare')
    parser_prepare.add_argument('--input', required=True)
    parser_prepare.add_argument('--root', default='qe_response')
    parser_prepare.add_argument('--dim', type=int, choices=[0, 3], default=3)
    parser_prepare.add_argument('--periodic-axes', default=None)
    parser_prepare.add_argument('--tr2-ph', type=float, default=1.0e-14)
    parser_prepare.add_argument('--no-raman', action='store_true')
    parser_prepare.add_argument('--force', action='store_true')

    parser_run = action_subparsers.add_parser('run')
    parser_run.add_argument('--root', default='qe_response')
    parser_run.add_argument('--pw-command', default='pw.x')
    parser_run.add_argument('--ph-command', default='ph.x')
    parser_run.add_argument('--dynmat-command', default='dynmat.x')
    parser_run.add_argument('--min-gap', type=float, default=0.01)
    parser_run.add_argument('--omp-threads', type=int, default=1)
    parser_run.add_argument('--dry-run', action='store_true')

    parser_status = action_subparsers.add_parser('status')
    parser_status.add_argument('--root', default='qe_response')

    parser_collect = action_subparsers.add_parser('collect')
    parser_collect.add_argument('--root', default='qe_response')
    parser_collect.add_argument('--broadening', type=float, default=8.0)
    parser_collect.add_argument('--points', type=int, default=2001)
    parser_collect.add_argument('--no-plot', action='store_true')

    parser_script = action_subparsers.add_parser('script')
    parser_script.add_argument('--root', default='qe_response')
    parser_script.add_argument(
        '--backend', choices=['shell', 'slurm', 'torque'], default='shell'
    )
    parser_script.add_argument('--output', default=None)
    parser_script.add_argument('--job-name', default='zstar-qe-response')
    parser_script.add_argument('--nodes', type=int, default=1)
    parser_script.add_argument('--tasks', type=int, default=1)
    parser_script.add_argument('--cpus-per-task', type=int, default=1)
    parser_script.add_argument('--walltime', default='24:00:00')
    parser_script.add_argument('--queue', default=None)
    parser_script.add_argument('--account', default=None)
    parser_script.add_argument('--env-script', default=None)


def _run_qe_action(args) -> None:
    from .qe_backend import (
        collect_qe_response,
        format_qe_status,
        generate_qe_backend_script,
        prepare_qe_response,
        qe_response_status,
        run_qe_response,
    )

    if args.qe_action == 'prepare':
        root = prepare_qe_response(
            args.input,
            args.root,
            dimensionality=args.dim,
            periodic_axes=args.periodic_axes,
            tr2_ph=args.tr2_ph,
            raman=not args.no_raman,
            force=args.force,
        )
        print(f"[OUT] {root}")
    elif args.qe_action == 'run':
        states = run_qe_response(
            args.root,
            pw_command=args.pw_command,
            ph_command=args.ph_command,
            dynmat_command=args.dynmat_command,
            min_gap_eV=args.min_gap,
            omp_threads=args.omp_threads,
            dry_run=args.dry_run,
        )
        print(format_qe_status(states))
        if any(state.status not in {'completed', 'dry-run'} for state in states):
            raise SystemExit(1)
    elif args.qe_action == 'status':
        print(format_qe_status(qe_response_status(args.root)))
    elif args.qe_action == 'collect':
        result = collect_qe_response(
            args.root,
            broadening_cm1=args.broadening,
            points=args.points,
            plot=not args.no_plot,
        )
        print(f"[OUT] {result['response_output']}")
    elif args.qe_action == 'script':
        output = generate_qe_backend_script(
            args.root,
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
        )
        print(f"[OUT] {output}")


def zstar_cli(argv=None) -> None:
    """
    Entry point function for the `zstar` command.

    If `argv` is None, arguments are taken from `sys.argv[1:]` (normal CLI use).
    """
    cli_argv = list(sys.argv[1:] if argv is None else argv)
    legacy_qe_command = None
    if cli_argv[:1] == ['qe']:
        legacy_qe_command = 'zstar qe'
        cli_argv[0] = 'qe-bec'
    elif cli_argv[:2] == ['backend', 'qe']:
        legacy_qe_command = 'zstar backend qe'
        cli_argv = ['qe-bec', *cli_argv[2:]]

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
    parser_gen.add_argument(
        '--dim', type=int, choices=[0, 1, 2, 3], default=3,
        help='Physical dimensionality: 0 for a molecule, 1 for a z-periodic wire, 2 for a slab, or 3 for bulk.'
    )
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
    gen_calc.add_argument('--cp2k', action='store_true',
                          help='Use a CP2K input template for Berry-phase BEC.')
    parser_gen.add_argument('--cp2k-root', default='cp2k_bec',
                            help='Output root for --cp2k (default: cp2k_bec).')
    parser_gen.add_argument('--displacement', type=float, default=0.01,
                            help='Atomic finite-displacement half-step in Angstrom.')

    # ---------------- deal ----------------
    parser_deal = subparsers.add_parser(
        'deal', help='Deal with polarization data to get BORN effective charge.'
    )
    parser_deal.add_argument(
        '--dim', type=int, choices=[0, 1, 2, 3], default=3,
        help='Physical dimensionality: 0 for a molecule, 1 for a z-periodic wire, 2 for a slab, or 3 for bulk.'
    )
    parser_deal.add_argument('--method', type=str, help='Finite difference method, by forward or central, with pricesion of first and seconde order. To save calculation resource you can choose forward', default='forward')
    parser_deal.add_argument(
        '--displacement', type=float, default=None,
        help='Finite-displacement half-step in Angstrom; default reads disp_Angstrom.out.'
    )
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
    deal_calc.add_argument('--cp2k', action='store_true',
                           help='Collect a generated CP2K Berry-phase BEC workflow.')
    parser_deal.add_argument('--cp2k-root', default='cp2k_bec')

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
    parser_born.add_argument('--dim', type=int, choices=[0, 1, 2, 3], default=3)
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
    parser_polar.add_argument('--dim', type=int, choices=[0, 1, 2, 3], default=3)
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
    polar_calc.add_argument('--cp2k', action='store_true',
                            help='Collect a generated CP2K Berry-phase BEC workflow.')
    parser_polar.add_argument('--cp2k-root', default='cp2k_bec')

    # ---------------- 2D polarization profile ----------------
    parser_polar2d = subparsers.add_parser(
        'polar2d',
        help='Visualize a cube-integrated out-of-plane slab polarization change.'
    )
    parser_polar2d.add_argument(
        '--reference-cube', required=True,
        help='Reference charge cube or directory containing one.'
    )
    parser_polar2d.add_argument(
        '--displaced-cube', required=True,
        help='Displaced charge cube or directory containing one.'
    )
    parser_polar2d.add_argument(
        '--displacement', type=float, default=None,
        help='Signed ionic displacement in Angstrom; enables an effective-charge report.'
    )
    parser_polar2d.add_argument(
        '--neutrality-tolerance', type=float, default=0.05
    )
    parser_polar2d.add_argument('--outdir', default='polarization_2d_profile')
    parser_polar2d.add_argument('--no-plot', action='store_true')

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
    parser_postph.add_argument(
        '--physical-dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_postph.add_argument(
        '--nac-model',
        choices=['bulk', 'gonze', 'wang', '1d-cutoff', '2d-cutoff'],
        default='gonze',
    )
    parser_postph.add_argument(
        '--q-direction', type=float, nargs=3, default=None,
        metavar=('QX', 'QY', 'QZ'),
    )
    parser_postph.add_argument('--symmprec', '--tol', type=float,
                               help='Symmetry precision of STRU, default is 1e-3',
                               default=1e-3)

    # ---------------- CP2K BEC backend ----------------
    parser_cp2k = subparsers.add_parser(
        'cp2k-bec', help='Prepare, run, and validate CP2K Born-charge workflows.'
    )
    cp2k_actions = parser_cp2k.add_subparsers(dest='cp2k_action', required=True)

    parser_cp2k_prepare = cp2k_actions.add_parser(
        'prepare', help='Generate reference and finite-displacement CP2K inputs.'
    )
    parser_cp2k_prepare.add_argument('--input', required=True)
    parser_cp2k_prepare.add_argument('--root', default='cp2k_bec')
    parser_cp2k_prepare.add_argument(
        '--dim', type=int, choices=[0, 3], default=3,
        help='0 for a molecular atomic polar tensor; 3 for periodic BEC.'
    )
    parser_cp2k_prepare.add_argument(
        '--method', choices=['forward', 'central'], default='central'
    )
    parser_cp2k_prepare.add_argument('--displacement', type=float, default=0.01)
    parser_cp2k_prepare.add_argument(
        '--atoms', default='all', help='One-based indices, e.g. 1,3-5, or all.'
    )
    parser_cp2k_prepare.add_argument('--force', action='store_true')

    parser_cp2k_run = cp2k_actions.add_parser(
        'run', help='Run CP2K displacement stages serially with restart support.'
    )
    parser_cp2k_run.add_argument('--root', default='cp2k_bec')
    parser_cp2k_run.add_argument('--cp2k-command', default='cp2k.psmp')
    parser_cp2k_run.add_argument('--omp-threads', type=int, default=1)
    parser_cp2k_run.add_argument('--data-dir', default=None)
    parser_cp2k_run.add_argument('--dry-run', action='store_true')
    parser_cp2k_run.add_argument('--stop-after', type=int, default=None)

    parser_cp2k_status = cp2k_actions.add_parser('status')
    parser_cp2k_status.add_argument('--root', default='cp2k_bec')

    parser_cp2k_collect = cp2k_actions.add_parser(
        'collect', help='Construct molecular APT or periodic BEC tensors from CP2K dipoles.'
    )
    parser_cp2k_collect.add_argument('--root', default='cp2k_bec')
    parser_cp2k_collect.add_argument('--output', default='Z-BORN-all.out')
    parser_cp2k_collect.add_argument('--json-output', default='cp2k_bec.json')
    parser_cp2k_collect.add_argument('--response-output', default='zstar_response.json')

    parser_cp2k_native = cp2k_actions.add_parser(
        'native', help='Run the CP2K 2025.2+ native finite-field APT reference.'
    )
    parser_cp2k_native.add_argument('--input', required=True)
    parser_cp2k_native.add_argument('--root', default='cp2k_native_apt')
    parser_cp2k_native.add_argument('--field-strength', type=float, default=3.0e-4)
    parser_cp2k_native.add_argument('--cp2k-command', default='cp2k.ssmp')
    parser_cp2k_native.add_argument('--omp-threads', type=int, default=1)
    parser_cp2k_native.add_argument('--data-dir', default=None)
    parser_cp2k_native.add_argument('--prepare-only', action='store_true')
    parser_cp2k_native.add_argument('--force', action='store_true')

    parser_cp2k_compare = cp2k_actions.add_parser(
        'compare', help='Compare ZStar dipole derivatives with CP2K native APT.'
    )
    parser_cp2k_compare.add_argument('--zstar-json', required=True)
    parser_cp2k_compare.add_argument('--native-apt', required=True)
    parser_cp2k_compare.add_argument('--output', default='cp2k_bec_comparison.json')

    # ---------------- VASP BEC backend ----------------
    parser_vasp_bec = subparsers.add_parser(
        'vasp-bec', help='Prepare, run, and collect VASP Born-charge workflows.'
    )
    vasp_bec_actions = parser_vasp_bec.add_subparsers(
        dest='vasp_bec_action', required=True
    )
    parser_vasp_bec_prepare = vasp_bec_actions.add_parser(
        'prepare', help='Prepare reference SCF and linear-response VASP stages.'
    )
    parser_vasp_bec_prepare.add_argument('--input-dir', default='.')
    parser_vasp_bec_prepare.add_argument('--root', default='vasp_bec')
    parser_vasp_bec_prepare.add_argument(
        '--method', choices=['dfpt', 'finite-field'], default='dfpt'
    )
    parser_vasp_bec_prepare.add_argument(
        '--field-strength', type=float, default=0.001,
        help='EFIELD_PEAD component in eV/Angstrom for finite-field mode.'
    )
    parser_vasp_bec_prepare.add_argument('--force', action='store_true')
    parser_vasp_bec_run = vasp_bec_actions.add_parser(
        'run', help='Run reference and response serially with an insulation gate.'
    )
    parser_vasp_bec_run.add_argument('--root', default='vasp_bec')
    parser_vasp_bec_run.add_argument('--vasp-command', default='vasp_std')
    parser_vasp_bec_run.add_argument('--min-gap', type=float, default=0.01)
    parser_vasp_bec_run.add_argument('--dry-run', action='store_true')
    parser_vasp_bec_status = vasp_bec_actions.add_parser('status')
    parser_vasp_bec_status.add_argument('--root', default='vasp_bec')
    parser_vasp_bec_collect = vasp_bec_actions.add_parser(
        'collect', help='Normalize OUTCAR BEC tensors and write ZStar/Phonopy files.'
    )
    parser_vasp_bec_collect.add_argument('--root', default='vasp_bec')
    parser_vasp_bec_collect.add_argument('--output', default='Z-BORN-all.out')
    parser_vasp_bec_collect.add_argument('--born-output', default='BORN')
    parser_vasp_bec_collect.add_argument('--json-output', default='vasp_bec.json')
    parser_vasp_bec_collect.add_argument('--response-output', default='zstar_response.json')
    parser_vasp_bec_script = vasp_bec_actions.add_parser(
        'script', help='Generate a shell, Slurm, or Torque serial driver.'
    )
    parser_vasp_bec_script.add_argument('--root', default='vasp_bec')
    parser_vasp_bec_script.add_argument(
        '--backend', choices=['shell', 'slurm', 'torque'], default='shell'
    )
    parser_vasp_bec_script.add_argument('--output', default=None)
    parser_vasp_bec_script.add_argument('--job-name', default='zstar-vasp-bec')
    parser_vasp_bec_script.add_argument('--nodes', type=int, default=1)
    parser_vasp_bec_script.add_argument('--tasks', type=int, default=1)
    parser_vasp_bec_script.add_argument('--cpus-per-task', type=int, default=1)
    parser_vasp_bec_script.add_argument('--walltime', default='24:00:00')
    parser_vasp_bec_script.add_argument('--queue', default=None)
    parser_vasp_bec_script.add_argument('--account', default=None)
    parser_vasp_bec_script.add_argument('--env-script', default=None)
    parser_vasp_bec_script.add_argument('--vasp-command', default=None)
    parser_vasp_bec_script.add_argument('--min-gap', type=float, default=0.01)
    parser_vasp_bec_compare = vasp_bec_actions.add_parser(
        'compare', help='Compare two normalized vasp_bec.json results.'
    )
    parser_vasp_bec_compare.add_argument('--first', required=True)
    parser_vasp_bec_compare.add_argument('--second', required=True)
    parser_vasp_bec_compare.add_argument('--output', default='vasp_bec_comparison.json')

    # ---------------- Quantum ESPRESSO BEC/response backend ----------------
    parser_qe = subparsers.add_parser(
        'qe-bec',
        help='Prepare, run, and collect Quantum ESPRESSO BEC/IR/Raman workflows.',
    )
    qe_actions = parser_qe.add_subparsers(dest='qe_action', required=True)
    _add_qe_action_parsers(qe_actions)

    # ---------------- calculator-native IR/Raman ----------------
    parser_spectra_backend = subparsers.add_parser(
        'spectra', help='Run calculator-native VASP or CP2K IR/Raman workflows.'
    )
    spectra_actions = parser_spectra_backend.add_subparsers(
        dest='spectra_action', required=True
    )
    parser_spectra_prepare = spectra_actions.add_parser(
        'prepare', help='Prepare a VASP or CP2K spectroscopy workflow.'
    )
    parser_spectra_prepare.add_argument(
        '--calculator', choices=['vasp', 'cp2k'], required=True
    )
    parser_spectra_prepare.add_argument('--root', default='calculator_spectra')
    parser_spectra_prepare.add_argument('--dim', type=int, choices=[0, 1, 2, 3], default=3)
    parser_spectra_prepare.add_argument('--input-dir', default='.')
    parser_spectra_prepare.add_argument('--modes-xml', default=None)
    parser_spectra_prepare.add_argument('--input', default=None)
    parser_spectra_prepare.add_argument('--modes', default=None)
    parser_spectra_prepare.add_argument('--acoustic-cutoff', type=float, default=5.0)
    parser_spectra_prepare.add_argument('--amplitude', type=float, default=0.02)
    parser_spectra_prepare.add_argument(
        '--method', choices=['dfpt', 'finite-field'], default='dfpt'
    )
    parser_spectra_prepare.add_argument('--field-strength', type=float, default=0.001)
    parser_spectra_prepare.add_argument('--cp2k-dx', type=float, default=0.01)
    parser_spectra_prepare.add_argument('--force', action='store_true')
    parser_spectra_run = spectra_actions.add_parser(
        'run', help='Run all prepared response stages serially and resumably.'
    )
    parser_spectra_run.add_argument('--root', default='calculator_spectra')
    parser_spectra_run.add_argument(
        '--command', dest='calculator_command', default=None,
        help='Calculator launch command (vasp_std or cp2k ... -o output.log).'
    )
    parser_spectra_run.add_argument('--min-gap', type=float, default=0.01)
    parser_spectra_run.add_argument('--omp-threads', type=int, default=1)
    parser_spectra_run.add_argument('--cp2k-data-dir', default=None)
    parser_spectra_run.add_argument('--dry-run', action='store_true')
    parser_spectra_run.add_argument('--stop-after', type=int, default=None)
    parser_spectra_status = spectra_actions.add_parser('status')
    parser_spectra_status.add_argument('--root', default='calculator_spectra')
    parser_spectra_script = spectra_actions.add_parser(
        'script', help='Generate a shell, Slurm, or Torque serial driver.'
    )
    parser_spectra_script.add_argument('--root', default='calculator_spectra')
    parser_spectra_script.add_argument(
        '--backend', choices=['shell', 'slurm', 'torque'], default='shell'
    )
    parser_spectra_script.add_argument('--output', default=None)
    parser_spectra_script.add_argument('--job-name', default='zstar-spectra')
    parser_spectra_script.add_argument('--nodes', type=int, default=1)
    parser_spectra_script.add_argument('--tasks', type=int, default=1)
    parser_spectra_script.add_argument('--cpus-per-task', type=int, default=1)
    parser_spectra_script.add_argument('--walltime', default='24:00:00')
    parser_spectra_script.add_argument('--queue', default=None)
    parser_spectra_script.add_argument('--account', default=None)
    parser_spectra_script.add_argument('--env-script', default=None)
    parser_spectra_script.add_argument('--command', dest='calculator_command', default=None)
    parser_spectra_script.add_argument('--min-gap', type=float, default=0.01)
    parser_spectra_collect = spectra_actions.add_parser(
        'collect', help='Collect both IR and Raman spectra and render plots.'
    )
    parser_spectra_collect.add_argument('--root', default='calculator_spectra')
    parser_spectra_collect.add_argument('--temperature', type=float, default=300.0)
    parser_spectra_collect.add_argument('--laser-nm', type=float, default=532.0)
    parser_spectra_collect.add_argument('--broadening', type=float, default=8.0)
    parser_spectra_collect.add_argument('--points', type=int, default=2001)
    parser_spectra_collect.add_argument('--no-plot', action='store_true')

    # ---------------- qNEP data bridge ----------------
    parser_qnep = subparsers.add_parser(
        'qnep', help='Prepare and audit GPUMD qNEP datasets with BEC labels.'
    )
    qnep_actions = parser_qnep.add_subparsers(dest='qnep_action', required=True)
    parser_qnep_augment = qnep_actions.add_parser(
        'augment', help='Append bec:R:9 labels to selected extxyz frames.'
    )
    parser_qnep_augment.add_argument('--input', required=True)
    parser_qnep_augment.add_argument('--output', default='train_qnep.xyz')
    qnep_source = parser_qnep_augment.add_mutually_exclusive_group(required=True)
    qnep_source.add_argument('--bec', help='BEC source for one frame.')
    qnep_source.add_argument('--map', dest='bec_map', help='CSV with frame,bec columns.')
    parser_qnep_augment.add_argument('--frame', type=int, default=0)
    parser_qnep_augment.add_argument('--audit-output', default=None)
    parser_qnep_check = qnep_actions.add_parser(
        'check', help='Validate qNEP extxyz structure and BEC columns.'
    )
    parser_qnep_check.add_argument('--input', required=True)
    parser_qnep_check.add_argument('--audit-output', default=None)
    parser_qnep_init = qnep_actions.add_parser(
        'init', help='Write a minimal qNEP nep.in from a labeled dataset.'
    )
    parser_qnep_init.add_argument('--input', required=True)
    parser_qnep_init.add_argument('--output', default='nep.in')
    parser_qnep_init.add_argument('--charge-mode', type=int, choices=[1, 2], default=2)
    parser_qnep_init.add_argument('--lambda-z', type=float, default=0.5)

    # ---------------- Born-charge database ----------------
    parser_db = subparsers.add_parser(
        'db', help='Collect ZStar workspaces into a High-K/BEC database.'
    )
    db_actions = parser_db.add_subparsers(dest='db_action', required=True)
    parser_db_init = db_actions.add_parser(
        'init', help='Write a candidate-manifest CSV template.'
    )
    parser_db_init.add_argument('--manifest', default='candidates.csv')
    parser_db_collect = db_actions.add_parser(
        'collect', help='Collect manifest workspaces into CSV and JSONL tables.'
    )
    parser_db_collect.add_argument('--manifest', default='candidates.csv')
    parser_db_collect.add_argument('--output', default='database')

    # ---------------- Agent Skill ----------------
    parser_agent_skill = subparsers.add_parser(
        'agent-skill',
        help='Install or inspect the packaged run-zstar-workflows Agent Skill.'
    )
    agent_skill_actions = parser_agent_skill.add_subparsers(
        dest='agent_skill_action', required=True
    )
    parser_agent_skill_install = agent_skill_actions.add_parser(
        'install', help='Install the packaged skill into an agent skills directory.'
    )
    parser_agent_skill_install.add_argument(
        '--dest', default=None,
        help='Parent skills directory; defaults to $CODEX_HOME/skills or ~/.codex/skills.'
    )
    parser_agent_skill_install.add_argument('--force', action='store_true')
    agent_skill_actions.add_parser('path', help='Print the packaged skill directory.')
    parser_agent_skill_preflight = agent_skill_actions.add_parser(
        'preflight', help='Emit a non-mutating JSON readiness report.'
    )
    parser_agent_skill_preflight.add_argument('--root', default='.')
    parser_agent_skill_preflight.add_argument(
        '--lane',
        choices=['bec', 'phonon', 'ir', 'raman', 'dielectric', 'md', 'cp2k', 'database'],
        default='bec',
    )
    parser_agent_skill_preflight.add_argument(
        '--dim', dest='dimensionality',
        choices=['molecule', '1d', '2d', 'bulk'], default='bulk'
    )

    # ---------------- calculator-neutral backend and response contracts ----------------
    parser_backend = subparsers.add_parser(
        'backend', help='Inspect calculator backend capabilities.'
    )
    backend_actions = parser_backend.add_subparsers(
        dest='backend_action', required=True
    )
    parser_backend_list = backend_actions.add_parser(
        'list', help='List implemented ZStar backend capabilities.'
    )
    parser_backend_list.add_argument('--json', action='store_true')
    parser_backend_list.add_argument(
        '--discover', action='store_true',
        help='Load third-party plugins registered under zstar.backends.'
    )
    parser_response = subparsers.add_parser(
        'response', help='Validate or convert calculator-neutral response records.'
    )
    response_actions = parser_response.add_subparsers(
        dest='response_action', required=True
    )
    parser_response_validate = response_actions.add_parser(
        'validate', help='Validate a zstar-response JSON document.'
    )
    parser_response_validate.add_argument('--input', required=True)
    parser_response_import = response_actions.add_parser(
        'import-bec', help='Convert an existing ZStar VASP/CP2K BEC JSON result.'
    )
    parser_response_import.add_argument('--input', required=True)
    parser_response_import.add_argument('--output', default='zstar_response.json')
    parser_response_import.add_argument(
        '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_response_import.add_argument(
        '--periodic-axes', default=None,
        help='Periodic axes, for example z, xy, or "x z".'
    )
    parser_response_abacus = response_actions.add_parser(
        'import-abacus', help='Normalize ABACUS/PYATB Z-BORN and optional BORN files.'
    )
    parser_response_abacus.add_argument('--zborn', required=True)
    parser_response_abacus.add_argument('--born', default=None)
    parser_response_abacus.add_argument('--output', default='zstar_response.json')
    parser_response_abacus.add_argument(
        '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_response_abacus.add_argument('--periodic-axes', default=None)
    parser_response_phonopy = response_actions.add_parser(
        'import-phonopy', help='Import Phonopy Gamma modes and optional BORN data.'
    )
    parser_response_phonopy.add_argument('--qpoints', default='qpoints.yaml')
    parser_response_phonopy.add_argument('--born', default=None)
    parser_response_phonopy.add_argument('--output', default='zstar_response.json')
    parser_response_phonopy.add_argument(
        '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_response_phonopy.add_argument('--periodic-axes', default=None)
    parser_response_intrinsic = response_actions.add_parser(
        'intrinsic', help='Add vacuum-independent 0D/1D/2D polarizability.'
    )
    parser_response_intrinsic.add_argument('--input', required=True)
    parser_response_intrinsic.add_argument('--output', default='zstar_response_intrinsic.json')
    parser_response_intrinsic.add_argument(
        '--lattice', type=float, nargs=9, required=True,
        metavar=('a1x', 'a1y', 'a1z', 'a2x', 'a2y', 'a2z', 'a3x', 'a3y', 'a3z'),
    )
    parser_response_intrinsic.add_argument(
        '--convention', choices=['gaussian', 'si-reduced'], default='gaussian'
    )

    parser_density = subparsers.add_parser(
        'density', help='Export calculator densities to the ZStar cube contract.'
    )
    density_actions = parser_density.add_subparsers(
        dest='density_action', required=True
    )
    parser_density_vasp = density_actions.add_parser('vasp-cube')
    parser_density_vasp.add_argument('--chgcar', required=True)
    parser_density_vasp.add_argument('--output', default='charge-density.cube')
    parser_density_vasp.add_argument('--potcar', default=None)
    parser_density_qe = density_actions.add_parser('qe-input')
    parser_density_qe.add_argument('--prefix', required=True)
    parser_density_qe.add_argument('--outdir', default='./tmp')
    parser_density_qe.add_argument('--cube', default='charge-density.cube')
    parser_density_qe.add_argument('--output', default='pp.in')
    parser_density_qe_sidecar = density_actions.add_parser('qe-sidecar')
    parser_density_qe_sidecar.add_argument('--cube', required=True)
    parser_density_qe_sidecar.add_argument('--pw-input', required=True)
    parser_density_qe_sidecar.add_argument('--pseudo-dir', required=True)
    parser_density_cp2k = density_actions.add_parser('cp2k-block')
    parser_density_cp2k.add_argument('--stride', type=int, nargs=3, default=(1, 1, 1))
    parser_density_cp2k.add_argument('--output', default='cp2k_density_cube.inc')
    parser_density_sidecar = density_actions.add_parser('sidecar')
    parser_density_sidecar.add_argument('--cube', required=True)
    parser_density_sidecar.add_argument('--backend', required=True)
    parser_density_sidecar.add_argument('--charges', type=float, nargs='+', required=True)

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
        '--dimensionality', '--dim', type=int, choices=[0, 1, 2, 3], default=3
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
        '--abacus-command', default=None,
        help='Override the backend-aware ABACUS launcher.'
    )
    parser_workflow_script.add_argument(
        '--pyatb-command', default=None,
        help='Override the backend-aware PYATB launcher.'
    )
    parser_workflow_script.add_argument('--mp-density', type=float, default=0.08)
    parser_workflow_script.add_argument(
        '--gap-mode', choices=['path', 'mp'], default='path'
    )
    parser_workflow_script.add_argument(
        '--dimensionality', '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_workflow_script.add_argument('--min-gap', type=float, default=0.01)
    parser_workflow_script.add_argument(
        '--no-insulation-check', action='store_true'
    )
    parser_workflow_script.add_argument(
        '--legacy-omega-max', type=float, default=30.0
    )
    parser_workflow_script.add_argument(
        '--dry-run', action='store_true',
        help='Generate a scheduler/environment smoke-test script without calculations.'
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
    parser_ir.add_argument('--dim', type=int, choices=[0, 1, 2, 3], default=3)
    parser_ir.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z',
        help='Periodic axis for --dim 1; default: z.'
    )
    parser_ir.add_argument(
        '--displacements', '--dipole-dir', default=None,
        help='For --dim 0, Raman-style +/- mode directory containing '
             'raman_manifest.json and PYATB polarizations.'
    )
    parser_ir.add_argument(
        '--polarization-subdir', default=None,
        help='Preferred PYATB polarization folder in each +/- stage '
             '(auto-detects pyatb and pyatb-polar by default).'
    )
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
    parser_raman_run.add_argument('--dim', type=int, choices=[0, 1, 2, 3], default=3)
    parser_raman_run.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z'
    )
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
    parser_raman_run.add_argument(
        '--ir-outdir', default='ir_spectrum',
        help='Molecular IR output directory used with --dim 0.'
    )
    parser_raman_run.add_argument('--no-spectrum', action='store_true')
    parser_raman_run.add_argument('--no-plot', action='store_true')
    parser_raman_run.add_argument(
        '--incident-polarization', type=float, nargs=3, default=None,
        metavar=('EX', 'EY', 'EZ'),
    )
    parser_raman_run.add_argument(
        '--scattered-polarization', type=float, nargs=3, default=None,
        metavar=('EX', 'EY', 'EZ'),
    )

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
        '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_raman_collect.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z'
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
        '--dim', type=int, choices=[0, 1, 2, 3], default=3
    )
    parser_raman_spectrum.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z'
    )
    parser_raman_spectrum.add_argument('--temperature', type=float, default=300.0)
    parser_raman_spectrum.add_argument('--laser-nm', type=float, default=532.0)
    parser_raman_spectrum.add_argument('--broadening', type=float, default=8.0)
    parser_raman_spectrum.add_argument('--max-frequency', type=float, default=None)
    parser_raman_spectrum.add_argument('--points', type=int, default=2001)
    parser_raman_spectrum.add_argument('--outdir', default='raman_spectrum')
    parser_raman_spectrum.add_argument('--no-plot', action='store_true')
    parser_raman_spectrum.add_argument(
        '--incident-polarization', type=float, nargs=3, default=None,
        metavar=('EX', 'EY', 'EZ'),
    )
    parser_raman_spectrum.add_argument(
        '--scattered-polarization', type=float, nargs=3, default=None,
        metavar=('EX', 'EY', 'EZ'),
    )

    parser_optics = subparsers.add_parser(
        'optics', help='Derive directional bulk optical constants from epsilon(omega).'
    )
    parser_optics.add_argument('--real', required=True)
    parser_optics.add_argument('--imag', required=True)
    parser_optics.add_argument(
        '--polarization', type=float, nargs=3, default=(1.0, 0.0, 0.0),
        metavar=('EX', 'EY', 'EZ'),
    )
    parser_optics.add_argument('--output', default='optical_constants.dat')

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
    parser_calc.add_argument('--dim', type=int, choices=[1, 2, 3], default=3)
    parser_calc.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z',
        help='Periodic axis for --dim 1; default: z.'
    )
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
    parser_freq.add_argument('--dim', type=int, choices=[1, 2, 3], default=3)
    parser_freq.add_argument(
        '--periodic-axis', choices=['x', 'y', 'z'], default='z',
        help='Periodic axis for --dim 1; default: z.'
    )
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
    md_bec.add_argument(
        '--bec-command',
        help='Batch predictor using ZSTAR_MD_REQUEST and ZSTAR_MD_OUTPUT NPZ/NPY contract.'
    )
    md_bec.add_argument(
        '--bec-provider',
        help='BEC provider entry-point name or module:object.'
    )
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
    parser_potential.add_argument('--vacuum-window', type=float, default=0.75,
                                  help='Local averaging width in Angstrom for each side-vacuum plateau.')
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

    args = parser.parse_args(cli_argv)

    if legacy_qe_command:
        print(
            f"[DEPRECATED] Use 'zstar qe-bec ...' instead of "
            f"'{legacy_qe_command} ...'.",
            file=sys.stderr,
        )

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
    if args.command == 'polar2d':
        from .polarization_2d import (
            compare_slab_charge_profiles,
            find_charge_cube,
            write_slab_charge_difference,
        )

        reference_cube = (
            find_charge_cube(args.reference_cube)
            if os.path.isdir(args.reference_cube)
            else args.reference_cube
        )
        displaced_cube = (
            find_charge_cube(args.displaced_cube)
            if os.path.isdir(args.displaced_cube)
            else args.displaced_cube
        )
        result = compare_slab_charge_profiles(
            reference_cube,
            displaced_cube,
            displacement_angstrom=args.displacement,
            neutrality_tolerance=args.neutrality_tolerance,
        )
        summary = write_slab_charge_difference(
            args.outdir, result, plot=not args.no_plot
        )
        print(
            "Slab dipole change: "
            f"{result.total_dipole_change_e_angstrom:.8g} e Angstrom"
        )
        if result.effective_charge_e is not None:
            print(f"Out-of-plane effective charge: {result.effective_charge_e:.8g} e")
        print(f"[OUT] {os.path.abspath(args.outdir)}")

    elif args.command == 'cp2k-bec':
        import json
        from pathlib import Path
        from .cp2k_bec import (
            collect_cp2k_bec,
            compare_cp2k_bec,
            cp2k_bec_status,
            format_cp2k_status,
            prepare_cp2k_bec,
            prepare_native_apt,
            run_cp2k_bec,
            run_native_apt,
        )

        if args.cp2k_action == 'prepare':
            root = prepare_cp2k_bec(
                args.input,
                args.root,
                method=args.method,
                displacement_angstrom=args.displacement,
                atoms=args.atoms,
                dimensionality=args.dim,
                force=args.force,
            )
            print(f"[OUT] {root}")
        elif args.cp2k_action == 'run':
            extra_env = {'CP2K_DATA_DIR': args.data_dir} if args.data_dir else None
            states = run_cp2k_bec(
                args.root,
                cp2k_command=args.cp2k_command,
                omp_threads=args.omp_threads,
                dry_run=args.dry_run,
                stop_after=args.stop_after,
                extra_env=extra_env,
            )
            print(format_cp2k_status(states))
        elif args.cp2k_action == 'status':
            print(format_cp2k_status(cp2k_bec_status(args.root)))
        elif args.cp2k_action == 'collect':
            result = collect_cp2k_bec(
                args.root,
                output=args.output,
                json_output=args.json_output,
                response_output=args.response_output,
            )
            print(f"[OUT] {result['output']}")
            print(f"[OUT] {result['json_output']}")
            print(f"[OUT] {result['response_output']}")
            sum_label = (
                "Acoustic-sum residual" if result["sum_scope"] == "all_atoms"
                else "Selected-atom tensor sum"
            )
            print(
                f"{sum_label} (e):\n"
                + "\n".join(
                    " ".join(f"{value: .6e}" for value in row)
                    for row in result['acoustic_sum_tensor']
                )
            )
        elif args.cp2k_action == 'native':
            root = prepare_native_apt(
                args.input,
                args.root,
                field_strength=args.field_strength,
                force=args.force,
            )
            print(f"[OUT] {root / 'input.inp'}")
            if not args.prepare_only:
                extra_env = {'CP2K_DATA_DIR': args.data_dir} if args.data_dir else None
                apt = run_native_apt(
                    root,
                    cp2k_command=args.cp2k_command,
                    omp_threads=args.omp_threads,
                    extra_env=extra_env,
                )
                print(f"[OUT] {apt}")
        elif args.cp2k_action == 'compare':
            result = compare_cp2k_bec(args.zstar_json, args.native_apt)
            output = Path(args.output).resolve()
            output.write_text(json.dumps(result, indent=2), encoding='utf-8')
            print(f"max |Delta tensor| = {result['max_abs']:.6e} e")
            print(f"RMS Delta tensor = {result['rms']:.6e} e")
            print(
                "max translational-sum residual: "
                f"ZStar={result['zstar_acoustic_sum_max_abs']:.6e} e, "
                f"native={result['native_acoustic_sum_max_abs']:.6e} e"
            )
            print(f"[OUT] {output}")

    elif args.command == 'vasp-bec':
        from .vasp_bec import (
            collect_vasp_bec,
            compare_vasp_bec,
            format_vasp_status,
            generate_vasp_backend_script,
            prepare_vasp_bec,
            run_vasp_bec,
            vasp_bec_status,
        )

        if args.vasp_bec_action == 'prepare':
            root = prepare_vasp_bec(
                args.input_dir,
                args.root,
                method=args.method,
                field_strength=args.field_strength,
                force=args.force,
            )
            print(f"[OUT] {root}")
        elif args.vasp_bec_action == 'run':
            states = run_vasp_bec(
                args.root,
                vasp_command=args.vasp_command,
                min_gap_eV=args.min_gap,
                dry_run=args.dry_run,
            )
            print(format_vasp_status(states))
        elif args.vasp_bec_action == 'status':
            print(format_vasp_status(vasp_bec_status(args.root)))
        elif args.vasp_bec_action == 'collect':
            result = collect_vasp_bec(
                args.root,
                output=args.output,
                born_output=args.born_output,
                json_output=args.json_output,
                response_output=args.response_output,
            )
            print(f"[OUT] {result['output']}")
            print(f"[OUT] {result['born_output']}")
            print(f"[OUT] {result['json_output']}")
            print(f"[OUT] {result['response_output']}")
            print(
                "Acoustic-sum residual (e):\n"
                + "\n".join(
                    " ".join(f"{value: .6e}" for value in row)
                    for row in result['acoustic_sum_tensor']
                )
            )
        elif args.vasp_bec_action == 'script':
            output = generate_vasp_backend_script(
                args.root,
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
                vasp_command=args.vasp_command,
                min_gap_eV=args.min_gap,
            )
            print(f"[OUT] {output}")
        elif args.vasp_bec_action == 'compare':
            import json
            from pathlib import Path

            result = compare_vasp_bec(args.first, args.second)
            output = Path(args.output).resolve()
            output.write_text(json.dumps(result, indent=2), encoding='utf-8')
            print(f"max |Delta Z*| = {result['bec_max_abs_e']:.6e} e")
            print(f"RMS Delta Z* = {result['bec_rms_e']:.6e} e")
            print(f"max |Delta epsilon_inf| = {result['epsilon_max_abs']:.6e}")
            print(f"[OUT] {output}")

    elif args.command == 'qe-bec':
        _run_qe_action(args)

    elif args.command == 'spectra':
        from .spectroscopy_backends import (
            calculator_spectra_status,
            collect_calculator_spectra,
            format_calculator_spectra_status,
            generate_calculator_spectra_script,
            prepare_cp2k_spectra,
            prepare_vasp_spectra,
            run_calculator_spectra,
        )

        if args.spectra_action == 'prepare':
            if args.calculator == 'vasp':
                if not args.modes_xml:
                    parser_spectra_prepare.error(
                        '--modes-xml is required for --calculator vasp'
                    )
                root = prepare_vasp_spectra(
                    args.input_dir,
                    args.modes_xml,
                    args.root,
                    amplitude=args.amplitude,
                    mode_numbers=_parse_modes(args.modes),
                    acoustic_cutoff_cm1=args.acoustic_cutoff,
                    method=args.method,
                    field_strength=args.field_strength,
                    dimensionality=args.dim,
                    force=args.force,
                )
            else:
                if not args.input:
                    parser_spectra_prepare.error(
                        '--input is required for --calculator cp2k'
                    )
                root = prepare_cp2k_spectra(
                    args.input,
                    args.root,
                    displacement_bohr=args.cp2k_dx,
                    dimensionality=args.dim,
                    force=args.force,
                )
            print(f"[OUT] {root}")
        elif args.spectra_action == 'run':
            extra_env = (
                {'CP2K_DATA_DIR': args.cp2k_data_dir}
                if args.cp2k_data_dir
                else None
            )
            states = run_calculator_spectra(
                args.root,
                command=args.calculator_command,
                min_gap_eV=args.min_gap,
                omp_threads=args.omp_threads,
                extra_env=extra_env,
                dry_run=args.dry_run,
                stop_after=args.stop_after,
            )
            print(format_calculator_spectra_status(states))
        elif args.spectra_action == 'status':
            print(
                format_calculator_spectra_status(
                    calculator_spectra_status(args.root)
                )
            )
        elif args.spectra_action == 'script':
            output = generate_calculator_spectra_script(
                args.root,
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
                calculator_command=args.calculator_command,
                min_gap_eV=args.min_gap,
            )
            print(f"[OUT] {output}")
        elif args.spectra_action == 'collect':
            result = collect_calculator_spectra(
                args.root,
                broadening_cm1=args.broadening,
                laser_nm=args.laser_nm,
                temperature_K=args.temperature,
                points=args.points,
                plot=not args.no_plot,
            )
            print(
                f"Collected {len(result['frequencies_cm-1'])} "
                f"{result['calculator'].upper()} modes."
            )
            print(f"[OUT] {os.path.abspath(args.root)}")

    elif args.command == 'qnep':
        from .qnep_dataset import (
            augment_qnep_dataset,
            check_qnep_dataset,
            write_qnep_input,
        )

        if args.qnep_action == 'augment':
            summary = augment_qnep_dataset(
                args.input,
                args.output,
                bec=args.bec,
                frame=args.frame,
                bec_map=args.bec_map,
                audit_output=args.audit_output,
            )
            print(
                f"Validated {summary['frames']} frames; "
                f"BEC labels added to {summary['labeled_frames']}."
            )
            print(f"[OUT] {summary['output']}")
            print(f"[OUT] {summary['audit_output']}")
        elif args.qnep_action == 'check':
            summary = check_qnep_dataset(
                args.input, audit_output=args.audit_output
            )
            print(
                f"Valid qNEP extxyz: frames={summary['frames']}, "
                f"labeled={summary['labeled_frames']}, "
                f"elements={','.join(summary['elements'])}"
            )
            if summary.get('audit_output'):
                print(f"[OUT] {summary['audit_output']}")
        elif args.qnep_action == 'init':
            output = write_qnep_input(
                args.input,
                args.output,
                charge_mode=args.charge_mode,
                lambda_z=args.lambda_z,
            )
            print(f"[OUT] {output}")

    elif args.command == 'db':
        from .bec_database import collect_database, write_manifest_template

        if args.db_action == 'init':
            output = write_manifest_template(args.manifest)
            print(f"[OUT] {output}")
        elif args.db_action == 'collect':
            summary = collect_database(args.manifest, args.output)
            print(
                "Collected "
                f"{summary['materials']} materials, {summary['atom_tensors']} atom tensors; "
                f"complete={summary['complete']}, incomplete={summary['incomplete']}, "
                f"rejected_metal={summary['rejected_metal']}."
            )
            print(f"[OUT] {os.path.abspath(args.output)}")

    elif args.command == 'agent-skill':
        from .agent_skill import (
            install_agent_skill,
            packaged_skill_path,
            write_preflight_json,
        )

        if args.agent_skill_action == 'install':
            destination = install_agent_skill(args.dest, force=args.force)
            print(f"Installed run-zstar-workflows to {destination}")
            print("Restart or open a new agent session to refresh skill discovery.")
        elif args.agent_skill_action == 'path':
            print(packaged_skill_path())
        elif args.agent_skill_action == 'preflight':
            print(
                write_preflight_json(
                    args.root,
                    lane=args.lane,
                    dimensionality=args.dimensionality,
                )
            )

    elif args.command == 'backend':
        import json

        from .backends import backend_capability_table, builtin_registry

        registry = builtin_registry()
        discovered = registry.discover() if args.discover else []
        if args.json:
            print(
                json.dumps(
                    {
                        'plugin_group': 'zstar.backends',
                        'discovered': discovered,
                        'backends': [backend.spec.to_dict() for backend in registry.list()],
                    },
                    indent=2,
                )
            )
        else:
            print(backend_capability_table(registry))
            if discovered:
                print(f"Discovered plugins: {', '.join(discovered)}")

    elif args.command == 'response':
        import json
        from pathlib import Path

        from .response_schema import (
            response_record_from_abacus_files,
            response_record_from_bec_result,
            validate_response_document,
        )

        if args.response_action == 'validate':
            print(json.dumps(validate_response_document(args.input), indent=2))
        elif args.response_action == 'import-bec':
            source = Path(args.input).resolve()
            data = json.loads(source.read_text(encoding='utf-8'))
            record = response_record_from_bec_result(
                data,
                dimensionality=args.dim,
                periodic_axes=args.periodic_axes,
                provenance={'source_file': str(source)},
            )
            output = record.write(args.output)
            print(f"[OUT] {output}")
        elif args.response_action == 'import-abacus':
            record = response_record_from_abacus_files(
                args.zborn,
                born_path=args.born,
                dimensionality=args.dim,
                periodic_axes=args.periodic_axes,
            )
            output = record.write(args.output)
            print(f"[OUT] {output}")
        elif args.response_action == 'import-phonopy':
            from .interoperability import response_record_from_phonopy

            record = response_record_from_phonopy(
                args.qpoints,
                born_path=args.born,
                dimensionality=args.dim,
                periodic_axes=args.periodic_axes,
            )
            output = record.write(args.output)
            print(f"[OUT] {output}")
        elif args.response_action == 'intrinsic':
            import numpy as np

            from .interoperability import add_intrinsic_response
            from .response_schema import ResponseRecord

            record = ResponseRecord.read(args.input)
            normalized = add_intrinsic_response(
                record,
                np.asarray(args.lattice, dtype=float).reshape(3, 3),
                convention=args.convention,
            )
            output = normalized.write(args.output)
            print(f"[OUT] {output}")

    elif args.command == 'density':
        from pathlib import Path

        from .density_adapters import (
            cp2k_density_cube_block,
            qe_pp_cube_input,
            vasp_chgcar_to_cube,
            write_cube_sidecar,
            write_qe_cube_sidecar,
        )

        if args.density_action == 'vasp-cube':
            output = vasp_chgcar_to_cube(
                args.chgcar, args.output, potcar_path=args.potcar
            )
        elif args.density_action == 'qe-input':
            output = Path(args.output).resolve()
            output.write_text(
                qe_pp_cube_input(
                    prefix=args.prefix,
                    outdir=args.outdir,
                    output_cube=args.cube,
                ),
                encoding='utf-8',
                newline='\n',
            )
        elif args.density_action == 'qe-sidecar':
            output = write_qe_cube_sidecar(
                args.cube, args.pw_input, pseudo_dir=args.pseudo_dir
            )
        elif args.density_action == 'cp2k-block':
            output = Path(args.output).resolve()
            output.write_text(
                cp2k_density_cube_block(stride=tuple(args.stride)),
                encoding='utf-8',
                newline='\n',
            )
        elif args.density_action == 'sidecar':
            output = write_cube_sidecar(
                args.cube, args.charges, backend=args.backend
            )
        print(f"[OUT] {output}")

    elif args.command == 'workflow':
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
                dry_run=args.dry_run,
            )
            print(f"[OUT] {script}")
            if args.submit:
                print(submit_backend_script(script, args.backend))

    elif args.command == 'ir':
        from .spectra import (
            calculate_molecular_ir_spectrum,
            calculate_ir_spectrum,
            collect_molecular_dipole_derivatives,
            load_gamma_modes,
            read_born_data,
            write_ir_outputs,
            write_molecular_ir_outputs,
        )

        modes = load_gamma_modes(args.qpoints)
        if args.dim == 0:
            if not args.displacements:
                parser_ir.error(
                    "--dim 0 requires --displacements/--dipole-dir"
                )
            mode_numbers, derivatives, _ = collect_molecular_dipole_derivatives(
                args.displacements,
                cell_volume_angstrom3=modes.volume_angstrom3,
                cell_lattice_angstrom=modes.lattice_angstrom,
                polarization_subdir=args.polarization_subdir,
            )
            result = calculate_molecular_ir_spectrum(
                modes,
                mode_numbers,
                derivatives,
                broadening_cm1=args.broadening,
                max_frequency_cm1=args.max_frequency,
                points=args.points,
            )
            summary = write_molecular_ir_outputs(
                args.outdir, result, plot=not args.no_plot
            )
        else:
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
                periodic_axis='xyz'.index(args.periodic_axis),
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
        import numpy as np

        from .spectra import (
            calculate_molecular_ir_spectrum,
            calculate_raman_spectrum,
            collect_molecular_dipole_derivatives,
            collect_raman_tensors,
            load_gamma_modes,
            load_raman_tensors,
            prepare_raman_displacements,
            write_raman_outputs,
            write_molecular_ir_outputs,
            write_native_line_spectrum_outputs,
        )
        from .spectroscopy_analysis import calculate_polarized_raman_spectrum

        incident_polarization = getattr(args, 'incident_polarization', None)
        scattered_polarization = getattr(args, 'scattered_polarization', None)
        if bool(incident_polarization) != bool(scattered_polarization):
            parser_raman.error(
                '--incident-polarization and --scattered-polarization are required together'
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
                molecular_ir=args.dim == 0,
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
                cell_volume_angstrom3=modes.volume_angstrom3,
                cell_cross_section_angstrom2=modes.cross_section_angstrom2(
                    'xyz'.index(args.periodic_axis)
                ),
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
                if args.incident_polarization:
                    polarized = calculate_polarized_raman_spectrum(
                        modes.frequencies_cm1[np.asarray(mode_numbers) - 1],
                        tensors,
                        mode_numbers=mode_numbers,
                        incident_polarization=args.incident_polarization,
                        scattered_polarization=args.scattered_polarization,
                        temperature_K=args.temperature,
                        laser_nm=args.laser_nm,
                        broadening_cm1=args.broadening,
                        max_frequency_cm1=args.max_frequency,
                        points=args.points,
                    )
                    polarized_dir = os.path.join(
                        args.spectrum_outdir, 'polarized'
                    )
                    write_native_line_spectrum_outputs(
                        polarized_dir,
                        polarized,
                        stem='raman_polarized',
                        plot=not args.no_plot,
                    )
                    print(f"[OUT] {os.path.abspath(polarized_dir)}")
                if args.dim == 0:
                    ir_numbers, derivatives, _ = (
                        collect_molecular_dipole_derivatives(
                            args.raman_dir,
                            cell_volume_angstrom3=modes.volume_angstrom3,
                            cell_lattice_angstrom=modes.lattice_angstrom,
                            polarization_subdir="pyatb-polar",
                        )
                    )
                    ir_result = calculate_molecular_ir_spectrum(
                        modes,
                        ir_numbers,
                        derivatives,
                        broadening_cm1=args.broadening,
                        max_frequency_cm1=args.max_frequency,
                        points=args.points,
                    )
                    ir_summary = write_molecular_ir_outputs(
                        args.ir_outdir,
                        ir_result,
                        plot=not args.no_plot,
                    )
                    print(f"Calculated {ir_summary['modes']} molecular IR modes.")
                    print(f"[OUT] {os.path.abspath(args.ir_outdir)}")
        elif args.raman_action == 'collect':
            mode_numbers, tensors, tensor_kind = collect_raman_tensors(
                args.raman_dir,
                dimensionality=args.dim,
                cell_height_angstrom=modes.cell_height_angstrom,
                cell_volume_angstrom3=modes.volume_angstrom3,
                cell_cross_section_angstrom2=modes.cross_section_angstrom2(
                    'xyz'.index(args.periodic_axis)
                ),
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
                    cell_volume_angstrom3=modes.volume_angstrom3,
                    cell_cross_section_angstrom2=modes.cross_section_angstrom2(
                        'xyz'.index(args.periodic_axis)
                    ),
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
            if args.incident_polarization:
                polarized = calculate_polarized_raman_spectrum(
                    modes.frequencies_cm1[np.asarray(mode_numbers) - 1],
                    tensors,
                    mode_numbers=mode_numbers,
                    incident_polarization=args.incident_polarization,
                    scattered_polarization=args.scattered_polarization,
                    temperature_K=args.temperature,
                    laser_nm=args.laser_nm,
                    broadening_cm1=args.broadening,
                    max_frequency_cm1=args.max_frequency,
                    points=args.points,
                )
                polarized_dir = os.path.join(args.outdir, 'polarized')
                write_native_line_spectrum_outputs(
                    polarized_dir,
                    polarized,
                    stem='raman_polarized',
                    plot=not args.no_plot,
                )
                print(f"[OUT] {os.path.abspath(polarized_dir)}")

    elif args.command == 'optics':
        from .spectroscopy_analysis import (
            optical_constants_from_dielectric,
            read_dielectric_response,
            write_optical_constants,
        )

        frequency, dielectric = read_dielectric_response(args.real, args.imag)
        result = optical_constants_from_dielectric(
            frequency, dielectric, polarization=args.polarization
        )
        output = write_optical_constants(args.output, result)
        print(f"[OUT] {output}")

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
            periodic_axis='xyz'.index(args.periodic_axis),
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
            periodic_axis='xyz'.index(args.periodic_axis),
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
            bec_command=args.bec_command,
            bec_provider=args.bec_provider,
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
            vacuum_window=args.vacuum_window,
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
        if getattr(args, 'cp2k', False):
            if not args.input:
                parser.error('zstar gen --cp2k requires --input CP2K_INPUT')
            from .cp2k_bec import prepare_cp2k_bec

            root = prepare_cp2k_bec(
                args.input,
                args.cp2k_root,
                method=args.method,
                displacement_angstrom=args.displacement,
                atoms=args.atom or 'all',
                dimensionality=args.dim,
                force=args.force,
            )
            print(f"[OUT] {root}")
            return

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
            method=method_fd,
            displacement_angstrom=args.displacement,
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
            dim=args.dim,
            physical_dim=args.physical_dim,
            nac_model=args.nac_model,
            q_direction=None if args.q_direction is None else tuple(args.q_direction),
        )

    elif args.command in ('rpolar', 'deal', 'born', 'polar'):
        if getattr(args, 'cp2k', False):
            from .cp2k_bec import collect_cp2k_bec

            result = collect_cp2k_bec(args.cp2k_root)
            print(f"[OUT] {result['output']}")
            print(f"[OUT] {result['json_output']}")
            return

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
            running_type=running_type,
            displacement_angstrom=getattr(args, 'displacement', None),
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
