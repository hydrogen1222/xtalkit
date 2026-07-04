"""CLI entry point for xtalkit."""

import argparse
import json
import os
import sys

from xtalkit import __version__
from xtalkit.spacegroup import (
    wyckoff_positions, sg_name, default_cell_params, crystal_system,
)
from xtalkit.marker import mark
from xtalkit.skeleton import generate
from xtalkit.enumerator import enumerate_structures
from xtalkit.builder import (
    AtomSite, FracAtom, find_wyckoff, free_params, build_structure,
    build_structure_frac, stoichiometry, stoichiometry_frac,
    format_formula, validate_cell, validate_atoms, validate_atoms_frac,
    detect_wyckoff,
)
from xtalkit.exporter import write_structure_cif, write_structure_xyz


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
    """Parse '4a:Xe,16e:Kr' into a dict (validated)."""
    from xtalkit.utils import parse_element_map
    return parse_element_map(value)


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
    # Verifies the bundled Wyckoff data is intact. Only 38 space groups are
    # populated so far (1-2, 195-230); the rest raise NotImplementedError,
    # which we count as "not yet supported" rather than a data error.
    try:
        supported = 0
        for n in range(1, 231):
            try:
                wyckoff_positions(n)
                supported += 1
            except NotImplementedError:
                continue
        print(f"[OK] Space group data intact ({supported}/230 space groups supported).")
    except Exception as e:
        print(f"[ERR] Space group data error: {e}")
        return 1
    return 0


def _parse_atom_raw(value: str) -> tuple[str, str, list[float]]:
    """Parse 'Li 16e 0.25' or 'Li 16e 0.25 1.0' into (element, wyckoff, numbers).

    The trailing numbers are free-coordinate values in template order, with an
    optional occupancy as the last value (resolved against the position's
    degree of freedom in ``cmd_build``).
    """
    parts = value.split()
    if len(parts) < 2:
        raise ValueError(
            f"--atom expects 'element wyckoff [free...] [occ]', got {value!r}"
        )
    element, wyckoff = parts[0], parts[1]
    try:
        numbers = [float(p) for p in parts[2:]]
    except ValueError as e:
        raise ValueError(f"--atom numeric parse error in {value!r}: {e}") from e
    return element, wyckoff, numbers


def _parse_atom_frac_raw(value: str) -> tuple[str, float, float, float, float]:
    """Parse 'Li 0.2563 0.2718 0.1832' or '... 0.6875' (with occupancy).

    Format: element x y z [occ]. Occupancy defaults to 1.0.
    """
    parts = value.split()
    if len(parts) < 4:
        raise ValueError(
            f"--atom-frac expects 'element x y z [occ]', got {value!r}"
        )
    element = parts[0]
    try:
        coords = [float(p) for p in parts[1:4]]
        occ = float(parts[4]) if len(parts) >= 5 else 1.0
    except ValueError as e:
        raise ValueError(f"--atom-frac numeric parse error in {value!r}: {e}") from e
    if not (0.0 <= occ <= 1.0 + 1e-6):
        raise ValueError(f"--atom-frac occupancy {occ} out of range [0, 1] in {value!r}")
    return element, coords[0], coords[1], coords[2], occ


def cmd_build(args) -> int:
    """Run the 'build' subcommand: assemble a CIF from refinement parameters."""
    try:
        if args.spec:
            with open(args.spec, encoding="utf-8") as f:
                spec = json.load(f)
            sg_number = int(spec["sg"])
            cell_params = spec["cell"]
            atoms = [
                AtomSite(
                    a["element"], a["wyckoff"],
                    a.get("free", []), float(a.get("occ", 1.0)),
                )
                for a in spec["atoms"]
            ]
            if not atoms:
                raise ValueError("No atoms specified.")
            for w in validate_cell(sg_number, cell_params) + \
                    validate_atoms(sg_number, atoms):
                print(f"[warn] {w}", file=sys.stderr)
            structure = build_structure(sg_number, cell_params, atoms)
            formula = format_formula(stoichiometry(sg_number, atoms))
        else:
            if args.sg is None:
                raise ValueError("--sg is required (or use --spec).")
            if args.cell is None:
                raise ValueError("--cell is required (or use --spec).")
            sg_number = args.sg
            cell_params = _parse_cell(args.cell)

            if args.atom_frac:
                # Fractional-coordinate mode: element x y z [occ] per atom.
                parsed = [_parse_atom_frac_raw(raw) for raw in args.atom_frac]
                # Drop occ==0 atoms (absent from the structure); the refinement
                # table may list them (e.g. a split site fully occupied by the
                # other species), so tolerate them with a note.
                frac_atoms = []
                for fa in parsed:
                    if fa[4] == 0.0:
                        print(f"[info] skipping {fa[0]} at "
                              f"({fa[1]},{fa[2]},{fa[3]}) — occupancy 0",
                              file=sys.stderr)
                    else:
                        frac_atoms.append(FracAtom(*fa))
                if not frac_atoms:
                    raise ValueError("at least one --atom-frac with occ > 0 is required.")
                for w in validate_cell(sg_number, cell_params) + \
                        validate_atoms_frac(sg_number, frac_atoms):
                    print(f"[warn] {w}", file=sys.stderr)
                # Report detected Wyckoff sites (helps catch input errors).
                print(f"[info] SG #{sg_number} ({sg_name(sg_number)}), "
                      f"{crystal_system(sg_number)} — detected Wyckoff sites:")
                for fa in frac_atoms:
                    label, mult = detect_wyckoff(sg_number, fa.x, fa.y, fa.z)
                    print(f"        {fa.element:<3} at ({fa.x:.4f},{fa.y:.4f},"
                          f"{fa.z:.4f}) occ={fa.occ:<5} -> {label} (mult {mult})")
                structure = build_structure_frac(sg_number, cell_params, frac_atoms)
                formula = format_formula(stoichiometry_frac(sg_number, frac_atoms))
            else:
                if not args.atom:
                    raise ValueError(
                        "at least one --atom or --atom-frac is required.")
                atoms = []
                for raw in args.atom:
                    element, wyckoff, numbers = _parse_atom_raw(raw)
                    wp = find_wyckoff(sg_number, wyckoff)
                    n_free = len(free_params(wp.coordinates))
                    if len(numbers) == n_free:
                        free, occ = numbers, 1.0
                    elif len(numbers) == n_free + 1:
                        free, occ = numbers[:-1], numbers[-1]
                    else:
                        raise ValueError(
                            f"--atom {raw!r}: {wyckoff} ({wp.coordinates!r}) "
                            f"expects {n_free} free value(s) plus optional occ, "
                            f"got {len(numbers)} number(s)."
                        )
                    atoms.append(AtomSite(element, wyckoff, free, occ))
                if not atoms:
                    raise ValueError("No atoms specified.")
                for w in validate_cell(sg_number, cell_params) + \
                        validate_atoms(sg_number, atoms):
                    print(f"[warn] {w}", file=sys.stderr)
                structure = build_structure(sg_number, cell_params, atoms)
                formula = format_formula(stoichiometry(sg_number, atoms))

        print(f"[OK] SG #{sg_number} ({sg_name(sg_number)}), "
              f"crystal system: {crystal_system(sg_number)}, formula: {formula}")
        output_base = args.output or f"SG{sg_number}_built"
        formats = [f.strip() for f in args.format.split(",") if f.strip()]
        paths = []
        for fmt in formats:
            path = f"{output_base}.{fmt}"
            if fmt == "cif":
                write_structure_cif(structure, path)
            elif fmt == "xyz":
                write_structure_xyz(structure, sg_number, path)
            else:
                raise ValueError(f"Unknown format {fmt!r} (use 'cif' and/or 'xyz').")
            paths.append(path)
        print(f"     Saved to: {', '.join(paths)}")
        return 0
    except (ValueError, FileNotFoundError, KeyError, NotImplementedError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_enumerate(args) -> int:
    """Run the 'enumerate' subcommand."""
    try:
        if args.min_cell_size > args.max_cell_size:
            raise ValueError(
                f"--min-cell-size ({args.min_cell_size}) must be <= "
                f"--max-cell-size ({args.max_cell_size})."
            )

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
            jobs=args.jobs,
            batch_size=args.batch_size,
            scratch_dir=args.scratch_dir,
            skip_preflight=args.skip_preflight,
        )
        print(f"[OK] Enumerated {len(paths)} structure(s).")
        print(f"     Output directory: {os.path.dirname(paths[0]) if paths else output_dir}")
        for p in paths[:10]:
            print(f"     - {p}")
        if len(paths) > 10:
            print(f"     ... and {len(paths) - 10} more")
        return 0
    except (RuntimeError, FileNotFoundError, ValueError) as e:
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
                        help="Matching tolerance in fractional coordinate units (default: 0.5)")
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
        help="Enumerate ordered configurations of a disordered CIF "
             "(requires the 'enumerate' extra + compiled enumlib)",
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
    p_enum.add_argument("--jobs", type=int, default=1,
                        help="Parallel workers for structure generation (default: 1 = "
                             "serial streaming; 0 = auto/cpu_count). Note: enum.x itself "
                             "is single-threaded and not sped up by this.")
    p_enum.add_argument("--batch-size", type=int, default=256,
                        help="Structures per batch — caps peak memory (default: 256). "
                             "Larger = fewer makestr.x calls but more memory.")
    p_enum.add_argument("--scratch-dir", type=str, default=None,
                        help="Directory for enumlib scratch files (default: system temp). "
                             "/dev/shm puts them on tmpfs (faster, limited to ~half RAM).")
    p_enum.add_argument("--skip-preflight", action="store_true", default=False,
                        help="Skip the non-integer-stoichiometry pre-check. DANGEROUS: "
                             "non-integer occupancy can make enumlib exhaust memory and "
                             "crash the system instead of failing cleanly.")
    p_enum.set_defaults(func=cmd_enumerate)

    # build
    p_build = sub.add_parser(
        "build",
        help="Build a CIF from refinement parameters (SG + cell + Wyckoff sites)",
    )
    p_build.add_argument("--sg", type=int, default=None,
                         help="Space group number (1-230)")
    p_build.add_argument("--cell", type=str, default=None,
                         help="Cell params: 'a b c alpha beta gamma'")
    p_build.add_argument("--atom", action="append", default=None, metavar="SPEC",
                         help="Atom (Wyckoff-letter mode): 'element wyckoff [free...] [occ]' "
                              "(repeatable), e.g. 'Li 16e 0.25' or 'Li 16e 0.25 1.0'")
    p_build.add_argument("--atom-frac", action="append", default=None, metavar="SPEC",
                         dest="atom_frac",
                         help="Atom (fractional-coord mode): 'element x y z [occ]' "
                              "(repeatable) — give refinement fractional coords directly, "
                              "e.g. 'Li 0.2563 0.2718 0.1832 0.6875'. Wyckoff orbit auto-detected.")
    p_build.add_argument("--spec", type=str, default=None,
                         help="JSON spec file (alternative to --sg/--cell/--atom)")
    p_build.add_argument("--format", type=str, default="cif",
                         help="Output format(s): cif,xyz (comma-separated, default: cif)")
    p_build.add_argument("-o", "--output", type=str, default=None,
                         help="Output base path (without extension)")
    p_build.set_defaults(func=cmd_build)

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
