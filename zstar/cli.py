# zstar/cli.py
# -*- coding: utf-8 -*-

"""
Unified CLI for the ZStar toolkit.

Subcommands:
- gen, deal, born, polar
- ph, postph
- wyckoff, irrep, vasp
- calc
- symcheck, bornsym

This is adapted from the original pykappa-dev CLI, with imports rewritten to
use the installed zstar package modules and exposed as the `zstar` entry point.
"""

import argparse
import os
import sys

from . import __version__
from .read_irrep import main as run_read_irrep_cli
from .calc_kappa import deal_q_vector as run_calc_kappa
from .gen_polar import gen_polar as run_gen
from .deal_polar import main as run_deal_polar
from .get_wyckoff import get_wyckoff_position as run_get_wyckoff
from .get_wyckoff import stru2vasp as run_stru2vasp
from .phonon_gen import run_phonopy_and_process_files as run_ph
from .phonon_post import run_eigen_irrep as run_postph
from .verify_born_symmetry import run_symcheck

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
    parser_gen.add_argument('--dim', type=int, help='Dim of your systems, 2 for 2D.', default=3)
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
                          help='Use PyATB for NSCF Berry phase.')

    # ---------------- deal ----------------
    parser_deal = subparsers.add_parser(
        'deal', help='Deal with polarization data to get BORN effective charge.'
    )
    parser_deal.add_argument('--dim', type=int,
                             help='Dim of your systems, 2 for 2D.', default=3)
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
                           help='Use PyATB for NSCF Berry phase.')

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
                           help='Use PyATB for NSCF Berry phase.')

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
                            help='Use PyATB for NSCF Berry phase.')

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

    p_ir = subparsers.add_parser('irrep', help='Classify Γ irreps from irreps.yaml')
    p_ir.add_argument('-f', '--file', default='irreps.yaml', help='Path to irreps.yaml')
    p_ir.add_argument('--mode', default='db', choices=['db', 'default', 'smodes'],
                      help='Classification mode: db/default (no external tools) or '
                           'smodes (legacy)')
    p_ir.add_argument('--acoustic-thz', type=float, default=0.05,
                      help='Acoustic threshold (THz)')
    p_ir.add_argument('--stru', default=None,
                      help='(compat) STRU path; smodes will look in CWD anyway')

    parser_vasp = subparsers.add_parser(
        'vasp', help='Convert ABACUS structure format STRU to VASP format POSCAR.'
    )
    parser_vasp.add_argument('--stru', help='Path to the STRU file', default='STRU')

    parser_calc = subparsers.add_parser('calc', help='Calculate dielectric tensor (kappa).')
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

    # ---------------- dispatch ----------------
    if args.command == 'irrep':
        run_read_irrep_cli(_build_irrep_argv(args))

    elif args.command == 'calc':
        run_calc_kappa(
            zero_tolerance=args.tolerance,
            ir_tolerance=args.ir_tolerance,
            ir_choose=args.ir_choose,
            plot_switch=args.plot,
            mode=args.mode,
            stru_file=args.stru_file,
            irreps_file=args.irreps_file
        )

    elif args.command == 'gen':
        # -------- gen 参数预处理：dim=2 时默认只沿 xy 方向位移 --------
        if args.move is None or str(args.move).strip() == "":
            if args.dim == 2:
                args.move = "x y"
                print("[INFO] dim=2 and --move not specified; default to move along 'x y' only.")
            else:
                args.move = "x y z"
                print("[INFO] dim=3 and --move not specified; default to move along 'x y z'.")

        move_input = [c for c in str(args.move) if c in ('x', 'y', 'z')]
        print("处理后的 --move 参数:", move_input)

        nscf_calculator = 'abacus' if getattr(args, 'abacus', False) else (
                          'pyatb' if getattr(args, 'pyatb', False) else 'pyatb')
        input_mode = args.input_mode or nscf_calculator
        input_sets = args.input_sets

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
            extract_starred_atoms_only=args.reduce
        )

    elif args.command == 'ph':
        run_ph(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            dim=args.dim,
            abacus_sub="abacus_x.sh",
            vasp_sub="vasp_scf.sh",
            node=args.node
        )

    elif args.command == 'postph':
        run_postph(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            nac=args.nac,
            dim=args.dim
        )

    elif args.command in ('rpolar', 'deal', 'born', 'polar'):
        calc_flag = 'abacus' if getattr(args, 'abacus', False) else 'pyatb'
        running_type = 'solo' if getattr(args, 'solo', False) else None

        kwargs = dict(
            f_stru=args.stru,
            symm_tol=args.symmprec,
            dimension=args.dim,
            running_type=running_type
        )
        if calc_flag:
            kwargs['nscf_calculator'] = calc_flag

        run_deal_polar(**kwargs)

    elif args.command == 'vasp':
        run_stru2vasp(f_stru=args.stru)

    elif args.command == 'wyckoff':
        run_get_wyckoff(fstru=args.stru)

    elif args.command == 'symcheck':
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
