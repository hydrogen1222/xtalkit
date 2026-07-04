"""Multiwfn-style TUI for xtalkit."""

import os
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from xtalkit.spacegroup import wyckoff_positions, sg_name, default_cell_params, crystal_system
from xtalkit.marker import mark
from xtalkit.skeleton import generate
from xtalkit.enumerator import enumerate_structures
from xtalkit.builder import (
    AtomSite, find_wyckoff, free_params, build_structure,
    stoichiometry, format_formula, validate_cell, validate_atoms,
)
from xtalkit.exporter import write_structure_cif, write_structure_xyz

console = Console()


def _header(text: str) -> None:
    """Print a styled header."""
    console.print(Panel(f"[bold cyan]{text}"), justify="center")


def _success(text: str) -> None:
    """Print a success message."""
    console.print(f"[green]✓ {text}[/green]")


def _error(text: str) -> None:
    """Print an error message."""
    console.print(f"[red]✗ {text}[/red]")


def _warn(text: str) -> None:
    """Print a warning."""
    console.print(f"[yellow]⚠ {text}[/yellow]")


def _prompt(text: str) -> str:
    """Prompt for user input."""
    return console.input(f"[bold]{text}[/bold] ").strip()


def resolve_cif_path(path: str) -> str | None:
    """Resolve a potentially relative CIF path. Returns abspath or None."""
    try:
        resolved = os.path.abspath(os.path.expanduser(path))
        if os.path.isfile(resolved):
            return resolved
        # Also try relative to cwd
        alt = os.path.join(os.getcwd(), path)
        if os.path.isfile(alt):
            return alt
        return None
    except Exception:
        return None


def select_wyckoff_positions(available: list[str]) -> list[str]:
    """Prompt user to select Wyckoff positions. Returns validated list."""
    while True:
        choice = _prompt("Wyckoff positions to mark (comma-separated, or 'all')")
        if choice.lower() == "all":
            return available
        parts = [p.strip() for p in choice.split(",")]
        invalid = [p for p in parts if p not in available]
        if invalid:
            _error(f"Invalid: {', '.join(invalid)}. Valid: {', '.join(available)}")
            continue
        return parts


def select_output_formats() -> list[str]:
    """Prompt for output format selection. Returns list of format strings."""
    fmt_map = {"1": ["cif"], "2": ["xyz"], "3": ["cif", "xyz"]}
    while True:
        choice = _prompt("Output format: [1] cif  [2] xyz  [3] all")
        if choice in fmt_map:
            return fmt_map[choice]
        _error("Invalid selection (1-3)")


def _prompt_yes_no(text: str, default: bool = False) -> bool:
    """Prompt for a yes/no choice."""
    suffix = "[Y/n]" if default else "[y/N]"
    while True:
        value = _prompt(f"{text} {suffix}")
        if not value:
            return default
        if value.lower() in {"y", "yes"}:
            return True
        if value.lower() in {"n", "no"}:
            return False
        _error("Please enter y or n")


def _prompt_int_or_none(text: str, default: int | None = None) -> int | None:
    """Prompt for an integer or None."""
    while True:
        value = _prompt(text)
        if not value:
            return default
        try:
            return int(value)
        except ValueError:
            _error("Please enter a number")


def _prompt_element_map() -> dict[str, str] | None:
    """Prompt for an element override map, validating format.

    Returns None if the user skips (empty input). Re-prompts on malformed
    input instead of crashing the workflow.
    """
    from xtalkit.utils import parse_element_map
    while True:
        map_str = _prompt("Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip")
        if not map_str:
            return None
        try:
            return parse_element_map(map_str)
        except ValueError as e:
            _error(str(e))


def _show_element_assignments(
    wyckoff_letters: list[str],
    element_map: dict[str, str] | None,
) -> None:
    """Display which dummy element is assigned to each Wyckoff letter."""
    from xtalkit.utils import assign_dummy_elements
    assignment = assign_dummy_elements(wyckoff_letters, element_map)
    console.print()
    table = Table(title="Dummy Element Assignments")
    table.add_column("Wyckoff", style="cyan")
    table.add_column("Element", style="yellow")
    for letter in sorted(
        wyckoff_letters,
        key=lambda w: (int("".join(c for c in w if c.isdigit()) or 0), w),
    ):
        table.add_row(letter, assignment[letter])
    console.print(table)


def _mark_workflow() -> None:
    """Interactive mark workflow."""
    _header("Mark Wyckoff Positions in CIF")

    # Get CIF path
    while True:
        path = _prompt("CIF file path")
        resolved = resolve_cif_path(path)
        if resolved:
            _success(f"found {resolved}")
            break
        _error(f"File not found: {path}")

    # Get SG number
    while True:
        try:
            sg_str = _prompt("Space group number")
            sg_number = int(sg_str)
            if 1 <= sg_number <= 230:
                break
            _error("Space group number must be 1-230")
        except ValueError:
            _error("Please enter a number")

    # Show available Wyckoff positions
    positions = wyckoff_positions(sg_number)
    console.print(f"\n  Space Group #{sg_number}: [cyan]{sg_name(sg_number)}[/cyan]")
    console.print(f"  Available Wyckoff positions: [yellow]{'  '.join(p.letter for p in positions)}[/yellow]\n")

    available = [p.letter for p in positions]
    selected = select_wyckoff_positions(available)

    # Mode
    while True:
        mode_choice = _prompt("Mode: [1] Overlay  [2] Replace")
        if mode_choice == "1":
            mode = "overlay"
            break
        elif mode_choice == "2":
            mode = "replace"
            break
        _error("Invalid selection (1-2)")

    # Element map
    element_map = _prompt_element_map()

    # Show dummy element assignments
    _show_element_assignments(selected, element_map)

    # Output format
    formats = select_output_formats()

    # Offset (default 0.02 for visibility in overlay, 0 for replace)
    default_offset = "0.02" if mode == "overlay" else "0"
    off_str = _prompt(f"Dummy atom offset in fractional coords (default {default_offset})")
    offset = float(off_str) if off_str else float(default_offset)

    # Tolerance
    tol_str = _prompt("Tolerance in fractional coords (default 0.5)")
    tolerance = float(tol_str) if tol_str else 0.5

    # Output base
    base = os.path.splitext(resolved)[0]
    out = _prompt(f"Output base path [default: {base}_WYCK]")
    output_base = out if out else f"{base}_WYCK"

    # Run
    try:
        result = mark(
            cif_path=resolved,
            sg_number=sg_number,
            wyckoff_letters=selected,
            mode=mode,
            tolerance=tolerance,
            element_map=element_map,
            formats=formats,
            output_base=output_base,
            offset=offset,
        )
        _success(f"Saved to: {result}")
    except Exception as e:
        _error(str(e))


def _skeleton_workflow() -> None:
    """Interactive skeleton generation workflow."""
    _header("Generate Wyckoff Skeleton")

    # Get SG number
    while True:
        try:
            sg_str = _prompt("Space group number")
            sg_number = int(sg_str)
            if 1 <= sg_number <= 230:
                break
            _error("Space group number must be 1-230")
        except ValueError:
            _error("Please enter a number")

    # Show Wyckoff positions
    positions = wyckoff_positions(sg_number)
    console.print(f"\n  Space Group #{sg_number}: [cyan]{sg_name(sg_number)}[/cyan]")
    console.print(f"  Available Wyckoff positions: [yellow]{'  '.join(p.letter for p in positions)}[/yellow]\n")

    available = [p.letter for p in positions]
    selected = select_wyckoff_positions(available)

    # Cell parameters
    defaults = default_cell_params(sg_number)
    console.print(
        f"\n  Default cell for SG #{sg_number}: "
        f"a={defaults['a']} b={defaults['b']} c={defaults['c']} "
        f"α={defaults['alpha']} β={defaults['beta']} γ={defaults['gamma']}"
    )
    while True:
        cell_choice = _prompt("[1] Use default  [2] Enter manually")
        if cell_choice == "1":
            cell_params = None
            break
        elif cell_choice == "2":
            cell_str = _prompt("Enter: a b c alpha beta gamma")
            parts = cell_str.split()
            if len(parts) == 6:
                cell_params = {
                    "a": float(parts[0]), "b": float(parts[1]), "c": float(parts[2]),
                    "alpha": float(parts[3]), "beta": float(parts[4]), "gamma": float(parts[5]),
                }
                break
            _error("Need 6 values: a b c alpha beta gamma")
        else:
            _error("Invalid selection (1-2)")

    # Output format
    formats = select_output_formats()

    # Output base
    out = _prompt(f"Output base path [default: SG{sg_number}_skeleton]")
    output_base = out if out else f"SG{sg_number}_skeleton"

    # Element map
    element_map = _prompt_element_map()

    # Run
    try:
        result = generate(
            sg_number=sg_number,
            wyckoff_letters=selected,
            cell_params=cell_params,
            element_map=element_map,
            formats=formats,
            output_base=output_base,
        )
        _success(f"Saved to: {result}")
    except Exception as e:
        _error(str(e))


def _info_workflow() -> None:
    """Interactive space group info workflow."""
    _header("Query Space Group Information")

    while True:
        try:
            sg_str = _prompt("Space group number")
            sg_number = int(sg_str)
            if 1 <= sg_number <= 230:
                break
            _error("Space group number must be 1-230")
        except ValueError:
            _error("Please enter a number")

    name = sg_name(sg_number)
    positions = wyckoff_positions(sg_number)
    cell = default_cell_params(sg_number)

    console.print(f"\n[bold]Space Group #{sg_number}: {name}[/bold]")
    console.print(f"Default cell: a={cell['a']} b={cell['b']} c={cell['c']} "
                  f"α={cell['alpha']} β={cell['beta']} γ={cell['gamma']}")

    table = Table(title=f"Wyckoff Positions ({len(positions)})")
    table.add_column("Letter", style="cyan")
    table.add_column("Multiplicity", style="yellow")
    table.add_column("Site Symmetry", style="green")
    table.add_column("Coordinates", style="white")

    for p in positions:
        table.add_row(p.letter, str(p.multiplicity), p.site_symmetry, p.coordinates)

    console.print(table)


def _fetch_workflow() -> None:
    """Interactive fetch / verify workflow."""
    from xtalkit.spacegroup import wyckoff_positions
    _header("Verify Space Group Database")
    try:
        supported = 0
        for n in range(1, 231):
            try:
                wyckoff_positions(n)
                supported += 1
            except NotImplementedError:
                continue
        _success(f"Space group data intact ({supported}/230 supported)")
    except Exception as e:
        _error(f"Data error: {e}")


def _build_workflow() -> None:
    """Interactive build workflow: assemble a CIF from refinement parameters."""
    _header("Build CIF from Refinement Parameters")

    # Space group
    while True:
        try:
            sg_str = _prompt("Space group number")
            sg_number = int(sg_str)
            if 1 <= sg_number <= 230:
                break
            _error("Space group number must be 1-230")
        except ValueError:
            _error("Please enter a number")

    # Wyckoff menu
    positions = wyckoff_positions(sg_number)
    console.print(f"\n  Space Group #{sg_number}: [cyan]{sg_name(sg_number)}[/cyan] "
                  f"({crystal_system(sg_number)})")
    table = Table(title=f"Wyckoff Positions ({len(positions)})")
    table.add_column("Label", style="cyan")
    table.add_column("Mult", style="yellow")
    table.add_column("Site Sym", style="green")
    table.add_column("Coordinates", style="white")
    for p in positions:
        table.add_row(p.letter, str(p.multiplicity), p.site_symmetry, p.coordinates)
    console.print(table)
    valid_labels = {p.letter for p in positions}

    # Cell parameters
    defaults = default_cell_params(sg_number)
    console.print(
        f"\n  Crystal system: {crystal_system(sg_number)}. Default cell: "
        f"a={defaults['a']} b={defaults['b']} c={defaults['c']} "
        f"α={defaults['alpha']} β={defaults['beta']} γ={defaults['gamma']}"
    )
    while True:
        cell_choice = _prompt("[1] Use default  [2] Enter manually")
        if cell_choice == "1":
            cell_params = defaults
            break
        elif cell_choice == "2":
            cell_str = _prompt("Enter: a b c alpha beta gamma")
            parts = cell_str.split()
            if len(parts) == 6:
                cell_params = dict(zip(
                    ["a", "b", "c", "alpha", "beta", "gamma"],
                    [float(x) for x in parts]))
                break
            _error("Need 6 values: a b c alpha beta gamma")
        else:
            _error("Invalid selection (1-2)")
    for w in validate_cell(sg_number, cell_params):
        _warn(w)

    # Add atoms
    atoms: list[AtomSite] = []
    while True:
        console.print()
        element = _prompt("Element symbol (e.g. Li, Na, Fe)")
        while True:
            label = _prompt("Wyckoff label (e.g. 4a, 16e)")
            if label in valid_labels:
                break
            _error(f"Invalid label. Valid: {', '.join(sorted(valid_labels))}")
        wp = find_wyckoff(sg_number, label)
        fp = free_params(wp.coordinates)
        free_vals: list[float] = []
        if fp:
            console.print(f"  {label} coordinates: [white]{wp.coordinates}[/white] "
                          f"— free parameters: {', '.join(fp)}")
            for var in fp:
                while True:
                    v = _prompt(f"  value for {var}")
                    try:
                        free_vals.append(float(v))
                        break
                    except ValueError:
                        _error("Please enter a number")
        else:
            console.print(f"  {label}: [white]{wp.coordinates}[/white] "
                          f"(no free parameters)")
        occ_str = _prompt("Occupancy (default 1.0)")
        try:
            occ = float(occ_str) if occ_str else 1.0
        except ValueError:
            _warn("Invalid occupancy, using 1.0")
            occ = 1.0
        atoms.append(AtomSite(element, label, free_vals, occ))
        _success(f"Added {element} on {label} (occ {occ}). Composition so far: "
                 f"{format_formula(stoichiometry(sg_number, atoms))}")
        for w in validate_atoms(sg_number, atoms):
            _warn(w)
        if _prompt("Add another atom? [y/N]").lower() != "y":
            break

    if not atoms:
        _error("No atoms added; nothing to build.")
        return

    # Output
    formats = select_output_formats()
    out = _prompt(f"Output base path [default: SG{sg_number}_built]")
    output_base = out if out else f"SG{sg_number}_built"

    try:
        structure = build_structure(sg_number, cell_params, atoms)
        _success(f"Formula: {format_formula(stoichiometry(sg_number, atoms))}")
        for fmt in formats:
            path = f"{output_base}.{fmt}"
            if fmt == "cif":
                write_structure_cif(structure, path)
            elif fmt == "xyz":
                write_structure_xyz(structure, sg_number, path)
            _success(f"Saved: {path}")
    except Exception as e:  # noqa: BLE001
        _error(str(e))


def _enumerate_workflow() -> None:
    """Interactive enumeration workflow."""
    _header("Enumerate Ordered Configurations")

    # Get CIF path
    while True:
        path = _prompt("CIF file path (with partial/disordered occupancy)")
        resolved = resolve_cif_path(path)
        if resolved:
            _success(f"found {resolved}")
            break
        _error(f"File not found: {path}")

    # Min cell size
    min_str = _prompt("Min supercell size (default 1)")
    try:
        min_cell = int(min_str) if min_str else 1
    except ValueError:
        _error("Invalid number, using default 1")
        min_cell = 1

    # Max cell size
    max_str = _prompt("Max supercell size (default 2)")
    try:
        max_cell = int(max_str) if max_str else 2
    except ValueError:
        _error("Invalid number, using default 2")
        max_cell = 2
    if max_cell < min_cell:
        _warn("Max < Min, adjusting max = min")
        max_cell = min_cell

    # Symmetry precision
    sp_str = _prompt("Symmetry tolerance (default 0.1)")
    try:
        symm_prec = float(sp_str) if sp_str else 0.1
    except ValueError:
        symm_prec = 0.1

    # Vacancy symbol
    vac = _prompt("Vacancy DummySpecies symbol (default X)") or "X"

    # Output dir
    base = os.path.splitext(os.path.basename(resolved))[0]
    out_str = _prompt(f"Output directory (default {base}_enum/)")
    output_dir = out_str if out_str else f"{base}_enum"

    # Max structures
    ms_str = _prompt("Max structures to write (Enter for unlimited)")
    try:
        max_structures = int(ms_str) if ms_str else None
    except ValueError:
        max_structures = None

    # Output format
    fmt_choice = _prompt("Output format: [1] cif  [2] xyz (default cif)")
    out_format = "xyz" if fmt_choice == "2" else "cif"

    # Timeout (none by default, but expose it so a runaway enumlib run can
    # be bounded interactively instead of hanging the TUI forever)
    to_str = _prompt("Timeout in minutes (Enter for none)")
    try:
        timeout = float(to_str) if to_str else None
    except ValueError:
        _warn("Invalid timeout, using none")
        timeout = None

    # Performance options (memory + parallelism for the structure-generation
    # phase; enum.x itself stays single-threaded regardless).
    console.print(
        "\n  [dim]Performance: --jobs parallelises structure generation "
        "(not enum.x); --batch-size caps memory; --scratch-dir /dev/shm "
        "avoids disk I/O.[/dim]"
    )
    jobs_str = _prompt("Parallel jobs (default 1 = serial; 0 = auto)")
    try:
        jobs = int(jobs_str) if jobs_str else 1
    except ValueError:
        _warn("Invalid jobs, using 1 (serial)")
        jobs = 1
    bs_str = _prompt("Batch size (default 256; smaller = less memory)")
    try:
        batch_size = int(bs_str) if bs_str else 256
    except ValueError:
        _warn("Invalid batch size, using 256")
        batch_size = 256
    scratch = _prompt("Scratch dir for enumlib (Enter for system temp; try /dev/shm)")
    scratch_dir = scratch if scratch else None

    console.print()
    console.print("[cyan]Running enumlib... (may take a few minutes)[/cyan]")
    try:
        paths = enumerate_structures(
            cif_path=resolved,
            min_cell_size=min_cell,
            max_cell_size=max_cell,
            symm_prec=symm_prec,
            vacancy_symbol=vac,
            output_dir=output_dir,
            max_structures=max_structures,
            timeout=timeout,
            format=out_format,
            jobs=jobs,
            batch_size=batch_size,
            scratch_dir=scratch_dir,
        )
        _success(f"Enumerated {len(paths)} structure(s)")

        # Try to show summary table
        try:
            import gemmi
            table = Table(title="Enumerated Structures")
            table.add_column("#", style="cyan")
            table.add_column("File", style="yellow")
            table.add_column("Formula", style="green")
            table.add_column("Atoms", style="white")
            table.add_column("a (Å)", style="white")
            for i, p in enumerate(paths):
                try:
                    doc = gemmi.cif.read_file(p)
                    block = doc.sole_block()
                    a = float(block.find_value("_cell_length_a") or 0)
                    atoms = list(block.find_values("_atom_site_label"))
                    formula = block.find_value("_chemical_formula_sum")
                    formula = str(formula).strip("'\"") if formula else "?"
                    table.add_row(str(i), os.path.basename(p), formula,
                                  str(len(atoms)), f"{a:.3f}")
                except Exception:
                    table.add_row(str(i), os.path.basename(p), "?", "?", "?")
            console.print(table)
        except ImportError:
            for i, p in enumerate(paths):
                console.print(f"  [{i}] {p}")

        _success(f"Saved to: {output_dir}")
    except Exception as e:
        _error(str(e))


def _ewald_workflow() -> None:
    """Interactive batch Ewald workflow."""
    _header("Batch Ewald Electrostatic Energy Ranking")

    while True:
        path = _prompt("CIF file, file folder, or root folder")
        resolved = resolve_cif_path(path)
        if resolved:
            _success(f"found {resolved}")
            break
        if os.path.isdir(os.path.abspath(os.path.expanduser(path))):
            resolved = os.path.abspath(os.path.expanduser(path))
            _success(f"found {resolved}")
            break
        _error(f"File or directory not found: {path}")

    while True:
        layout_choice = _prompt("Directory layout: [1] flat files  [2] one-level nested")
        if layout_choice == "1":
            layout = "flat"
            break
        if layout_choice == "2":
            layout = "nested"
            break
        _error("Invalid selection (1-2)")

    while True:
        energy_choice = _prompt(
            "Oxidation source: [1] explicit charges  [2] guess  [3] already present"
        )
        if energy_choice == "1":
            from xtalkit.ewald import parse_charges
            while True:
                raw = _prompt("Charges, e.g. Li:1 Ge:4 P:5 S:-2")
                try:
                    charges = parse_charges(raw.split())
                    guess = False
                    break
                except ValueError as e:
                    _error(str(e))
            break
        if energy_choice == "2":
            charges = None
            guess = True
            break
        if energy_choice == "3":
            charges = None
            guess = False
            break
        _error("Invalid selection (1-3)")

    per_atom = _prompt_yes_no("Rank by energy per atom", default=False)

    while True:
        sort_choice = _prompt("Sort order: [1] lowest first  [2] highest first")
        if sort_choice == "1":
            descending = False
            break
        if sort_choice == "2":
            descending = True
            break
        _error("Invalid selection (1-2)")

    top_n = _prompt_int_or_none("Keep top N structures after sorting (Enter for all)")
    if top_n is not None and top_n <= 0:
        _error("Top N must be a positive number")
        top_n = None

    group = _prompt_yes_no("Copy selected structures into folders", default=False)
    move = False
    if group:
        move = _prompt_yes_no("Move instead of copy", default=False)

    base = os.path.splitext(os.path.basename(resolved.rstrip(os.sep)))[0] or "ewald"
    out_dir = _prompt(f"Output directory [default: {base}_ewald]")
    output_dir = out_dir if out_dir else f"{base}_ewald"

    console.print()
    console.print("[cyan]Ranking structures...[/cyan]")
    try:
        from xtalkit.ewald import batch_ewald, group_rows, split_rows, write_csv

        rows = batch_ewald(
            [resolved],
            charges=charges,
            guess=guess,
            per_atom=per_atom,
            layout=layout,
        )
        rows = sorted(rows, key=lambda row: row.ewald_energy, reverse=descending)
        selected, remaining = split_rows(rows, top_n)

        table = Table(title=f"Ewald Ranking ({len(rows)} structures)")
        table.add_column("Rank", style="cyan")
        table.add_column("File", style="yellow")
        table.add_column("Formula", style="green")
        table.add_column("Atoms", style="white")
        table.add_column("Energy (eV)", style="white")
        preview = selected if top_n is not None else rows
        for i, row in enumerate(preview[:20], 1):
            table.add_row(
                str(i),
                row.relative_path,
                row.formula,
                str(row.n_atoms),
                f"{row.ewald_energy:.6f}",
            )
        console.print(table)
        if len(preview) > 20:
            console.print(f"[dim]... and {len(preview) - 20} more[/dim]")

        csv_path = os.path.join(output_dir, "ranking.csv")
        write_csv(rows, csv_path)
        _success(f"Saved ranking CSV: {csv_path}")

        if group:
            group_info = group_rows(
                selected,
                remaining,
                output_dir,
                move=move,
            )
            _success(
                f"Grouped into {group_info['selected_dir']} "
                f"and {group_info['remaining_dir']}"
            )
    except Exception as e:
        _error(str(e))


def _show_main_menu() -> None:
    """Display the main menu."""
    console.clear()
    console.print()
    panel = Panel(
        "[bold white]  [1] Mark CIF    -- Mark Wyckoff in a structure[/bold white]\n"
        "  [2] Skeleton    -- Generate pure Wyckoff skeleton\n"
        "  [3] Query SG    -- View space group information\n"
        "  [4] Fetch DB    -- Verify database online\n"
        "  [5] Enumerate   -- Enumerate ordered configurations (enumlib)\n"
        "  [6] Build       -- Build CIF from refinement params (SG+cell+sites)\n"
        "  [7] Ewald       -- Batch Ewald ranking and grouping\n"
        "  [0] Exit",
        title="xtalkit . Crystal Wyckoff Toolkit",
        border_style="cyan",
    )
    console.print(panel)
    console.print()


def run_tui() -> int:
    """Launch the interactive TUI."""
    workflows = {
        "1": _mark_workflow,
        "2": _skeleton_workflow,
        "3": _info_workflow,
        "4": _fetch_workflow,
        "5": _enumerate_workflow,
        "6": _build_workflow,
        "7": _ewald_workflow,
    }

    while True:
        _show_main_menu()
        choice = _prompt("Input your choice")

        if choice == "0":
            console.print("[cyan]Goodbye![/cyan]")
            return 0

        func = workflows.get(choice)
        if func:
            console.print()
            func()
            console.print()
            _prompt("Press Enter to continue...")
        else:
            _error("Invalid choice (0-7)")
