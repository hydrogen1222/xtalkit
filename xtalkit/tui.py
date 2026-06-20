"""Multiwfn-style TUI for xtalkit."""

import os
import sys

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint

from xtalkit.spacegroup import wyckoff_positions, sg_name, default_cell_params
from xtalkit.marker import mark
from xtalkit.skeleton import generate

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
    fmt_map = {"1": ["cif"], "2": ["vesta"], "3": ["xyz"], "4": ["cif", "vesta", "xyz"]}
    while True:
        choice = _prompt("Output format: [1] cif  [2] vesta  [3] xyz  [4] all")
        if choice in fmt_map:
            return fmt_map[choice]
        _error("Invalid selection (1-4)")


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

    # Output format
    formats = select_output_formats()

    # Tolerance
    tol_str = _prompt("Tolerance in A (default 0.5)")
    tolerance = float(tol_str) if tol_str else 0.5

    # Element map
    map_str = _prompt("Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip")
    element_map = None
    if map_str:
        element_map = {}
        for pair in map_str.split(","):
            letter, elem = pair.split(":")
            element_map[letter.strip()] = elem.strip()

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
        f"a={defaults['alpha']} b={defaults['beta']} g={defaults['gamma']}"
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
    map_str = _prompt("Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip")
    element_map = None
    if map_str:
        element_map = {}
        for pair in map_str.split(","):
            letter, elem = pair.split(":")
            element_map[letter.strip()] = elem.strip()

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
                  f"a={cell['alpha']} b={cell['beta']} g={cell['gamma']}")

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
    _header("Verify Space Group Database")
    try:
        for n in range(1, 231):
            wyckoff_positions(n)
        _success("Space group data intact (230/230 OK)")
    except Exception as e:
        _error(f"Data error: {e}")


def _show_main_menu() -> None:
    """Display the main menu."""
    console.clear()
    console.print()
    panel = Panel(
        "[bold white]  [1] Mark CIF    -- Mark Wyckoff in a structure[/bold white]\n"
        "  [2] Skeleton    -- Generate pure Wyckoff skeleton\n"
        "  [3] Query SG    -- View space group information\n"
        "  [4] Fetch DB    -- Verify database online\n"
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
            _error("Invalid choice (0-4)")
