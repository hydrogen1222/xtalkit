"""CLI entry point for xtalkit."""

import argparse
import os
import sys

from xtalkit import __version__
from xtalkit.spacegroup import (
    wyckoff_positions, sg_name, default_cell_params, crystal_system,
)
from xtalkit.marker import mark
from xtalkit.skeleton import generate
from xtalkit.enumerator import enumerate_structures


def _parse_wyckoff(value: str) -> list[str]:
    """Parse comma-separated Wyckoff letters."""
    parts = [p.strip() for p in value.split(",")]
    if not parts:
        raise argparse.ArgumentTypeError("At least one Wyckoff letter required")
    return parts


def _parse_cell(value: str) -> dict:
    """Parse 'a b c alpha beta gamma' into a dict."""
    parts = value.split()
    if len(parts) != 6:
        raise argparse.ArgumentTypeError(
            "Cell parameters must be: a b c alpha beta gamma"
        )
    keys = ["a", "b", "c", "alpha", "beta", "gamma"]
    return {k: float(v) for k, v in zip(keys, parts)}


def _parse_elems(value: str) -> dict[str, str]:
    """Parse '4a:Xe,16e:Kr' into a dict."""
    mapping = {}
    for pair in value.split(","):
        letter, elem = pair.split(":")
        mapping[letter.strip()] = elem.strip()
    return mapping


def cmd_mark(args) -> int:
    """Run the 'mark' subcommand."""
    try:
        # Expand "all" to all Wyckoff letters for this SG
        if args.wyckoff == ["all"]:
            from xtalkit.spacegroup import wyckoff_positions
            args.wyckoff = [w.letter for w in wyckoff_positions(args.sg)]

        # Build output base name
        if args.output:
            output_base = args.output
        else:
            base = os.path.splitext(os.path.basename(args.cif))[0]
            output_base = os.path.join(os.path.dirname(args.cif) or ".", f"{base}_WYCK")

        formats = args.format.split(",")
        element_map = _parse_elems(args.map) if args.map else None

        result = mark(
            cif_path=args.cif,
            sg_number=args.sg,
            wyckoff_letters=args.wyckoff,
            mode=args.mode,
            tolerance=args.tol,
            element_map=element_map,
            formats=formats,
            output_base=output_base,
            offset=args.offset,
        )
        print(f"[OK] Done. Saved to: {result}")
        return 0
    except (ValueError, FileNotFoundError, NotImplementedError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_skeleton(args) -> int:
    """Run the 'skeleton' subcommand."""
    try:
        # Expand "all" to all Wyckoff letters for this SG
        if args.wyckoff == ["all"]:
            from xtalkit.spacegroup import wyckoff_positions
            args.wyckoff = [w.letter for w in wyckoff_positions(args.sg)]

        if args.output:
            output_base = args.output
        else:
            output_base = f"SG{args.sg}_skeleton"

        formats = args.format.split(",")
        element_map = _parse_elems(args.map) if args.map else None
        cell_params = _parse_cell(args.cell) if args.cell else None

        result = generate(
            sg_number=args.sg,
            wyckoff_letters=args.wyckoff,
            cell_params=cell_params,
            element_map=element_map,
            formats=formats,
            output_base=output_base,
        )
        print(f"[OK] Done. Saved to: {result}")
        return 0
    except (ValueError, FileNotFoundError, NotImplementedError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_info(args) -> int:
    """Run the 'info' subcommand."""
    try:
        sg_number = args.sg
        name = sg_name(sg_number)
        positions = wyckoff_positions(sg_number)
        cell = default_cell_params(sg_number)

        print(f"Space Group #{sg_number}: {name}")
        print(f"Crystal System: {crystal_system(sg_number)}")
        print(f"Default cell: a={cell['a']} b={cell['b']} c={cell['c']} "
              f"α={cell['alpha']} β={cell['beta']} γ={cell['gamma']}")
        print(f"\nWyckoff Positions ({len(positions)}):")
        print(f"  {'Letter':<8} {'Mult':<6} {'Site Sym':<10} {'Coordinates'}")
        print(f"  {'-'*8} {'-'*6} {'-'*10} {'-'*20}")
        for p in positions:
            print(f"  {p.letter:<8} {p.multiplicity:<6} {p.site_symmetry:<10} {p.coordinates}")

        return 0
    except (ValueError, FileNotFoundError, NotImplementedError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_fetch(args) -> int:
    """Run the 'fetch' subcommand."""
    from xtalkit.spacegroup import wyckoff_positions
    # Currently a no-op: Gemmi data is bundled with the library.
    # Future: could check for Gemmi updates or pull from Bilbao.
    try:
        # Validate we can still query all 230 SGs
        for n in range(1, 231):
            wyckoff_positions(n)
        print("[OK] Space group data is intact (230/230 space groups OK).")
    except Exception as e:
        print(f"[ERR] Space group data error: {e}")
        return 1
    return 0


def cmd_enumerate(args) -> int:
    """Run the 'enumerate' subcommand."""
    try:
        base = os.path.splitext(os.path.basename(args.cif))[0]
        output_dir = args.output_dir or f"{base}_enum"

        paths = enumerate_structures(
            cif_path=args.cif,
            min_cell_size=args.min_cell_size,
            max_cell_size=args.max_cell_size,
            symm_prec=args.symm_prec,
            vacancy_symbol=args.vacancy_symbol,
            output_dir=output_dir,
            max_structures=args.max_structures,
            timeout=args.timeout,
            format=args.format,
        )
        print(f"[OK] Enumerated {len(paths)} structure(s).")
        print(f"     Output directory: {os.path.dirname(paths[0]) if paths else output_dir}")
        for p in paths[:10]:
            print(f"     - {p}")
        if len(paths) > 10:
            print(f"     ... and {len(paths) - 10} more")
        return 0
    except (RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Build the full CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="xtalkit",
        description="Crystal Wyckoff Toolkit — mark Wyckoff positions for VESTA visualization",
    )
    parser.add_argument("--version", action="version", version=f"xtalkit {__version__}")
    sub = parser.add_subparsers(dest="command")

    # mark
    p_mark = sub.add_parser("mark", help="Mark Wyckoff positions in a CIF file")
    p_mark.add_argument("cif", help="Path to the input CIF file")
    p_mark.add_argument("--sg", type=int, required=True, help="Space group number (1-230)")
    p_mark.add_argument("--wyckoff", type=_parse_wyckoff, required=True,
                        help="Wyckoff letters to mark (comma-separated, or 'all')")
    p_mark.add_argument("--mode", choices=["overlay", "replace"], default="overlay",
                        help="overlay: add dummies; replace: swap real atoms (default: overlay)")
    p_mark.add_argument("--tol", type=float, default=0.5,
                        help="Matching tolerance in angstrom (default: 0.5)")
    p_mark.add_argument("--offset", type=float, default=0.02,
                        help="Fractional offset for dummy atoms (default: 0.02, use 0 for exact positions)")
    p_mark.add_argument("--map", type=str, default=None,
                        help="Element override, e.g. '4a:Xe,16e:Kr'")
    p_mark.add_argument("--format", type=str, default="cif",
                        help="Output format(s): cif,xyz (comma-separated, default: cif)")
    p_mark.add_argument("-o", "--output", type=str, default=None,
                        help="Output base path (without extension)")
    p_mark.set_defaults(func=cmd_mark)

    # skeleton
    p_skel = sub.add_parser("skeleton", help="Generate pure Wyckoff skeleton")
    p_skel.add_argument("--sg", type=int, required=True, help="Space group number (1-230)")
    p_skel.add_argument("--wyckoff", type=_parse_wyckoff, required=True,
                        help="Wyckoff letters to generate (comma-separated, or 'all')")
    p_skel.add_argument("--cell", type=str, default=None,
                        help="Cell params: 'a b c alpha beta gamma'")
    p_skel.add_argument("--map", type=str, default=None,
                        help="Element override, e.g. '4a:Xe,16e:Kr'")
    p_skel.add_argument("--format", type=str, default="cif",
                        help="Output format(s): cif,xyz (comma-separated)")
    p_skel.add_argument("-o", "--output", type=str, default=None,
                        help="Output base path (without extension)")
    p_skel.set_defaults(func=cmd_skeleton)

    # info
    p_info = sub.add_parser("info", help="Query space group Wyckoff information")
    p_info.add_argument("--sg", type=int, required=True, help="Space group number (1-230)")
    p_info.set_defaults(func=cmd_info)

    # fetch
    p_fetch = sub.add_parser("fetch", help="Verify space group database")
    p_fetch.set_defaults(func=cmd_fetch)

    # enumerate
    p_enum = sub.add_parser(
        "enumerate",
        help="Enumerate ordered configurations of a disordered CIF (requires pymatgen + enumlib)",
    )
    p_enum.add_argument("cif", help="Path to the input disordered CIF")
    p_enum.add_argument("--min-cell-size", type=int, default=1,
                        help="Minimum supercell size to enumerate (default: 1)")
    p_enum.add_argument("--max-cell-size", type=int, default=2,
                        help="Maximum supercell size to enumerate (default: 2)")
    p_enum.add_argument("--symm-prec", type=float, default=0.1,
                        help="Symmetry tolerance for SpacegroupAnalyzer (default: 0.1)")
    p_enum.add_argument("--vacancy-symbol", type=str, default="X",
                        help="DummySpecies symbol for vacancies (default: X)")
    p_enum.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: <cif_basename>_enum/)")
    p_enum.add_argument("--max-structures", type=int, default=None,
                        help="Stop after N structures (default: unlimited)")
    p_enum.add_argument("--timeout", type=float, default=None,
                        help="Timeout in minutes (default: none)")
    p_enum.add_argument("--format", type=str, default="cif", choices=["cif", "xyz"],
                        help="Output format (default: cif)")
    p_enum.set_defaults(func=cmd_enumerate)

    return parser


def main() -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        # No subcommand: launch TUI
        from xtalkit.tui import run_tui
        return run_tui()
    else:
        return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
