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


def _parse_scaling_matrix(values: list[str] | None) -> list[list[int]]:
    """Parse 3 diagonal values or a full 3x3 scaling matrix."""
    if not values:
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    nums = [int(v) for v in values]
    if len(nums) == 3:
        return [[nums[0], 0, 0], [0, nums[1], 0], [0, 0, nums[2]]]
    if len(nums) == 9:
        return [nums[0:3], nums[3:6], nums[6:9]]
    raise ValueError("--scaling-matrix expects 3 diagonal values or 9 matrix values")


def _parse_charge_map(values: list[str] | None) -> dict[str, float]:
    """Parse species charge tokens like Li:1 S:-2."""
    charges = {}
    for token in values or []:
        if ":" not in token:
            raise ValueError(f"charge token {token!r} must be ELEMENT:CHARGE")
        element, charge = token.split(":", 1)
        charges[element] = float(charge)
    return charges


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


def cmd_shry_prepare(args) -> int:
    """Run 'shry prepare'."""
    try:
        from xtalkit.enumeration import prepare_shry_input
        result = prepare_shry_input(
            input_cif=args.cif,
            output_cif=args.out,
            vacancy_symbol=args.vacancy_symbol,
            occupancy_overrides=args.set_occupancy,
            parent_spacegroup=args.parent_spacegroup,
            target_formula=args.target_formula,
            strict=args.strict,
            symmetrize=args.symmetrize,
            symprec=args.symprec,
            angle_tolerance=args.angle_tolerance,
            scaling_matrix=_parse_scaling_matrix(args.scaling_matrix),
        )
        print(f"[OK] SHRY-ready CIF: {result['output_cif']}")
        print(f"     Manifest: {result['manifest']}")
        print(f"     Orbit grouping: {result['orbit_grouping']}")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_shry_count(args) -> int:
    """Run 'shry count'."""
    try:
        from xtalkit.enumeration import count_shry_structures
        result = count_shry_structures(
            cif_path=args.cif,
            scaling_matrix=_parse_scaling_matrix(args.scaling_matrix),
            symprec=args.symprec,
            angle_tolerance=args.angle_tolerance,
            atol=args.atol,
            strict=args.strict,
            symmetrize=args.symmetrize,
            out_json=args.out,
        )
        print(f"[OK] SHRY count-only: {result['count_only_result']}")
        print(f"     Saved to: {args.out}")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_shry_enum(args) -> int:
    """Run 'shry enum'."""
    try:
        from xtalkit.enumeration import enumerate_with_shry
        result = enumerate_with_shry(
            cif_path=args.cif,
            output_dir=args.out,
            scaling_matrix=_parse_scaling_matrix(args.scaling_matrix),
            symprec=args.symprec,
            angle_tolerance=args.angle_tolerance,
            atol=args.atol,
            expect_count=args.expect_count,
            remove_vacancy=args.remove_vacancy,
            target_formula=args.target_formula,
            write_cif_output=args.write_cif,
            write_poscar=args.write_poscar,
            write_degeneracy=args.write_degeneracy,
            dir_size=args.dir_size,
            strict=args.strict,
            symmetrize=args.symmetrize,
        )
        print(f"[OK] SHRY generated {result['generated_clean_count']} structure(s).")
        print(f"     Output directory: {args.out}")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_shry_verify(args) -> int:
    """Run 'shry verify'."""
    try:
        from xtalkit.enumeration import verify_shry_outputs
        result = verify_shry_outputs(
            output_dir=args.output_dir,
            check_count=args.check_count,
            check_formula=args.check_formula,
            check_dedup=args.check_dedup,
            target_formula=args.target_formula,
            vacancy_symbol=args.vacancy_symbol,
            symprec_list=[float(v) for v in args.symprec_list] if args.symprec_list else None,
            cross_backend=args.cross_backend,
            check_degeneracy=args.check_degeneracy,
        )
        print(f"[OK] SHRY verify passed for {result['clean_count']} structure(s).")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_shry_postprocess(args) -> int:
    """Run 'shry postprocess'."""
    try:
        from xtalkit.enumeration import (
            rank_by_shortest_distance,
            rank_by_ewald,
            write_tblite_inputs,
            write_cp2k_inputs,
            write_slurm_array,
        )
        outputs = []
        if args.shortest_distance:
            species = tuple(args.pair) if args.pair else None
            rows = rank_by_shortest_distance(args.output_dir, species=species)
            outputs.append(f"shortest-distance rows: {len(rows)}")
        if args.ewald:
            charges = _parse_charge_map(args.ewald_charges)
            rows = rank_by_ewald(args.output_dir, charges)
            outputs.append(f"ewald rows: {len(rows)}")
        if args.write_tblite:
            jobs = write_tblite_inputs(args.output_dir, charge=args.charge, uhf=args.uhf)
            outputs.append(f"tblite jobs: {len(jobs)}")
        if args.write_cp2k:
            jobs = write_cp2k_inputs(args.output_dir, template=args.cp2k_template)
            outputs.append(f"cp2k jobs: {len(jobs)}")
        if args.write_slurm:
            script = write_slurm_array(args.output_dir, job_kind=args.slurm_kind)
            outputs.append(f"slurm script: {script}")
        if not outputs:
            raise ValueError("no postprocess action selected")
        print("[OK] SHRY postprocess complete.")
        for line in outputs:
            print(f"     {line}")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def cmd_ewald(args) -> int:
    """Run 'ewald': batch Ewald electrostatic energy ranking."""
    try:
        from xtalkit.ewald import (
            batch_ewald,
            parse_charges,
            split_rows,
            group_rows,
            write_csv,
        )
        charges = parse_charges(args.charges) if args.charges else None
        base = os.path.splitext(os.path.basename(args.paths[0].rstrip(os.sep)))[0] or "ewald"
        group_dir = args.group_dir or f"{base}_ewald"
        out_csv = args.out or os.path.join(group_dir, "ranking.csv")
        rows = batch_ewald(
            args.paths,
            charges=charges,
            guess=args.guess,
            per_atom=args.per_atom,
            layout=args.layout,
        )
        descending = args.sort == "desc"
        rows = sorted(rows, key=lambda row: row.ewald_energy, reverse=descending)
        selected, remaining = split_rows(rows, args.top_n)
        energy_col = "E_per_atom (eV)" if args.per_atom else "Ewald E (eV)"
        order_label = "highest" if descending else "lowest"
        print(f"Ranked {len(rows)} structure(s) by Ewald energy "
              f"({order_label} first):")
        print(f"  {'Rank':<5} {'File':<32} {'Formula':<14} {'N':<5} {energy_col}")
        print(f"  {'-'*5} {'-'*32} {'-'*14} {'-'*5} {'-'*16}")
        preview_rows = selected if args.top_n is not None else rows
        for i, row in enumerate(preview_rows[:20], 1):
            print(f"  {i:<5} {row.relative_path[:32]:<32} {row.formula:<14} "
                  f"{row.n_atoms:<5} {row.ewald_energy:.6f}")
        if len(preview_rows) > 20:
            print(f"  ... and {len(preview_rows) - 20} more")
        write_csv(rows, out_csv)
        print(f"\nWrote ranking CSV: {out_csv}")
        if args.group:
            group_info = group_rows(
                selected,
                remaining,
                group_dir,
                selected_name=args.selected_name,
                remaining_name=args.remaining_name,
                move=args.move,
            )
            print(
                f"Grouped into: {group_info['selected_dir']} "
                f"and {group_info['remaining_dir']}"
            )
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Build the full CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="xtalkit",
        description="Crystal Wyckoff Toolkit — mark Wyckoff positions, build CIFs, and batch Ewald rank structures",
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

    # shry
    p_shry = sub.add_parser(
        "shry",
        help="Strict SHRY workflow for partially occupied structure enumeration",
    )
    shry_sub = p_shry.add_subparsers(dest="shry_command", required=True)

    p_shry_prepare = shry_sub.add_parser("prepare", help="Prepare a SHRY-ready CIF")
    p_shry_prepare.add_argument("cif", help="Input CIF with partial occupancy")
    p_shry_prepare.add_argument("--out", required=True, help="Output SHRY-ready CIF")
    p_shry_prepare.add_argument("--vacancy-symbol", default="X",
                                help="Pseudo species for vacancies (default: X)")
    p_shry_prepare.add_argument("--set-occupancy", default=None,
                                help="Explicit idealization, e.g. 'Li1:11/16 M1:Ge1/2P1/2'")
    p_shry_prepare.add_argument("--parent-spacegroup", type=int, default=None,
                                help="Expected parent space-group number")
    p_shry_prepare.add_argument("--target-formula", default=None,
                                help="Expected formula after removing vacancies")
    p_shry_prepare.add_argument("--scaling-matrix", nargs="+", default=None,
                                help="3 diagonal values or 9 full matrix values (default: 1 1 1)")
    p_shry_prepare.add_argument("--symprec", type=float, default=0.01)
    p_shry_prepare.add_argument("--angle-tolerance", type=float, default=5.0)
    p_shry_prepare.add_argument("--symmetrize", action="store_true", default=False)
    p_shry_prepare.add_argument(
        "--strict", action=argparse.BooleanOptionalAction, default=True,
        help="Strict audited workflow (default on). Use --no-strict to relax "
             "count/audit gating for exploratory runs.")
    p_shry_prepare.set_defaults(func=cmd_shry_prepare)

    p_shry_count = shry_sub.add_parser("count", help="Run SHRY count-only")
    p_shry_count.add_argument("cif", help="SHRY-ready CIF")
    p_shry_count.add_argument("--scaling-matrix", nargs="+", default=None,
                              help="3 diagonal values or 9 full matrix values (default: 1 1 1)")
    p_shry_count.add_argument("--symprec", type=float, default=0.01)
    p_shry_count.add_argument("--angle-tolerance", type=float, default=5.0)
    p_shry_count.add_argument("--atol", type=float, default=1e-5)
    p_shry_count.add_argument("--symmetrize", action="store_true", default=False)
    p_shry_count.add_argument("--out", default="count.json")
    p_shry_count.add_argument(
        "--strict", action=argparse.BooleanOptionalAction, default=True,
        help="Strict audited workflow (default on). Use --no-strict to relax.")
    p_shry_count.set_defaults(func=cmd_shry_count)

    p_shry_enum = shry_sub.add_parser("enum", help="Run strict SHRY enumeration")
    p_shry_enum.add_argument("cif", help="SHRY-ready CIF")
    p_shry_enum.add_argument("--scaling-matrix", nargs="+", default=None,
                             help="3 diagonal values or 9 full matrix values (default: 1 1 1)")
    p_shry_enum.add_argument("--symprec", type=float, default=0.01)
    p_shry_enum.add_argument("--angle-tolerance", type=float, default=5.0)
    p_shry_enum.add_argument("--atol", type=float, default=1e-5)
    p_shry_enum.add_argument("--expect-count", type=int, required=True)
    p_shry_enum.add_argument("--out", required=True, help="Output workflow directory")
    p_shry_enum.add_argument("--remove-vacancy", default="X")
    p_shry_enum.add_argument("--target-formula", default=None)
    p_shry_enum.add_argument("--write-cif", action="store_true", default=False)
    p_shry_enum.add_argument("--write-poscar", action="store_true", default=False)
    p_shry_enum.add_argument("--write-degeneracy", action="store_true", default=False)
    p_shry_enum.add_argument("--dir-size", type=int, default=10000)
    p_shry_enum.add_argument("--symmetrize", action="store_true", default=False)
    p_shry_enum.add_argument(
        "--strict", action=argparse.BooleanOptionalAction, default=True,
        help="Strict audited workflow (default on). Use --no-strict to relax "
             "the expect-count/mod-only gating.")
    p_shry_enum.set_defaults(func=cmd_shry_enum)

    p_shry_verify = shry_sub.add_parser("verify", help="Verify SHRY outputs")
    p_shry_verify.add_argument("output_dir")
    p_shry_verify.add_argument("--check-count", action="store_true", default=False)
    p_shry_verify.add_argument("--check-formula", action="store_true", default=False)
    p_shry_verify.add_argument("--check-dedup", action="store_true", default=False)
    p_shry_verify.add_argument("--target-formula", default=None)
    p_shry_verify.add_argument("--vacancy-symbol", default="X")
    p_shry_verify.add_argument("--symprec-list", nargs="*", default=None)
    p_shry_verify.add_argument("--cross-backend", choices=["supercell"], default=None)
    p_shry_verify.add_argument("--check-degeneracy", action="store_true", default=False)
    p_shry_verify.set_defaults(func=cmd_shry_verify)

    p_shry_post = shry_sub.add_parser("postprocess", help="Post-process SHRY outputs")
    p_shry_post.add_argument("output_dir")
    p_shry_post.add_argument("--shortest-distance", action="store_true", default=False)
    p_shry_post.add_argument("--pair", nargs=2, default=None,
                             help="Species pair for shortest-distance ranking, e.g. Li Li")
    p_shry_post.add_argument("--ewald", action="store_true", default=False)
    p_shry_post.add_argument("--ewald-charges", nargs="*", default=None,
                             help="Explicit charges, e.g. Li:1 Ge:4 P:5 S:-2")
    p_shry_post.add_argument("--write-tblite", action="store_true", default=False)
    p_shry_post.add_argument("--charge", type=int, default=0)
    p_shry_post.add_argument("--uhf", type=int, default=0)
    p_shry_post.add_argument("--write-cp2k", action="store_true", default=False)
    p_shry_post.add_argument("--cp2k-template", default=None)
    p_shry_post.add_argument("--write-slurm", action="store_true", default=False)
    p_shry_post.add_argument("--slurm-kind", choices=["tblite", "cp2k"], default="tblite")
    p_shry_post.set_defaults(func=cmd_shry_postprocess)

    # ewald
    p_ewald = sub.add_parser(
        "ewald",
        help="Batch Ewald electrostatic energy ranking for CIF/POSCAR folders",
    )
    p_ewald.add_argument("paths", nargs="+",
                         help="Input file(s) or directory root(s)")
    p_ewald.add_argument("--layout", choices=["flat", "nested"], default="flat",
                         help="Directory layout: flat files or one-level nested subdirs")
    p_ewald.add_argument("--sort", choices=["asc", "desc"], default="asc",
                         help="Sort order: asc = lowest first, desc = highest first")
    p_ewald.add_argument("--top-n", type=int, default=None,
                         help="Keep only the first N structures after sorting")
    p_ewald.add_argument("--charges", nargs="*", default=None,
                         help="Explicit oxidation states, e.g. Li:1 Ge:4 P:5 S:-2")
    p_ewald.add_argument("--guess", action="store_true", default=False,
                         help="Let pymatgen guess oxidation states")
    p_ewald.add_argument("--per-atom", action="store_true", default=False,
                         help="Rank by energy per atom instead of total energy")
    p_ewald.add_argument("--out", default=None,
                         help="Write full ranking CSV to this path")
    p_ewald.add_argument("--group", action="store_true", default=False,
                         help="Copy ranked structures into selected/rest folders")
    p_ewald.add_argument("--group-dir", default=None,
                         help="Grouping root directory (default: <input>_ewald)")
    p_ewald.add_argument("--selected-name", default="selected",
                         help="Folder name for the selected top-N structures")
    p_ewald.add_argument("--remaining-name", default="rest",
                         help="Folder name for the remaining structures")
    p_ewald.add_argument("--move", action="store_true", default=False,
                         help="Move files instead of copying them when grouping")
    p_ewald.set_defaults(func=cmd_ewald)

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
