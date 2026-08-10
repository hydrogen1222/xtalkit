"""Generate pure Wyckoff skeleton structures (no real atoms)."""

import gemmi
from xtalkit.spacegroup import wyckoff_positions, default_cell_params
from xtalkit.exporter import DummyAtom, write_cif, write_xyz
from xtalkit.builder import validate_cell
from xtalkit.utils import (
    assign_dummy_elements,
    representative_coord,
    unit_cell_from_parameters,
)


def generate(
    sg_number: int,
    wyckoff_letters: list[str],
    cell_params: dict | None,
    element_map: dict[str, str] | None,
    formats: list[str],
    output_base: str,
) -> str:
    """Generate a skeleton structure with only dummy atoms at Wyckoff positions.

    Returns comma-separated output file paths.
    """
    if not formats:
        raise ValueError("at least one output format is required")
    if len(set(formats)) != len(formats):
        raise ValueError("output formats must not contain duplicates")
    unsupported_formats = sorted(set(formats) - {"cif", "xyz"})
    if unsupported_formats:
        raise ValueError(
            f"unsupported output format(s): {', '.join(unsupported_formats)}"
        )
    if not wyckoff_letters:
        raise ValueError("at least one Wyckoff position is required")
    if len(set(wyckoff_letters)) != len(wyckoff_letters):
        raise ValueError("requested Wyckoff positions must not contain duplicates")

    # Get Wyckoff positions
    all_wyckoffs = wyckoff_positions(sg_number)
    valid_letters = {w.letter for w in all_wyckoffs}

    # Validate
    for letter in wyckoff_letters:
        if letter not in valid_letters:
            raise ValueError(
                f"Invalid Wyckoff letter '{letter}' for SG #{sg_number}. "
                f"Valid: {', '.join(sorted(valid_letters))}"
            )

    selected = [w for w in all_wyckoffs if w.letter in wyckoff_letters]
    assignment = assign_dummy_elements(wyckoff_letters, element_map)

    # Cell parameters
    if cell_params is None:
        params = default_cell_params(sg_number)
    else:
        params = cell_params
    cell = unit_cell_from_parameters(params)
    metric_issues = validate_cell(sg_number, params)
    if metric_issues:
        raise ValueError(
            "unit-cell metric is inconsistent with the selected space group: "
            + "; ".join(metric_issues)
        )

    # Build structure
    structure = gemmi.Structure()
    structure.cell = cell
    sg = gemmi.SpaceGroup(sg_number)
    structure.spacegroup_hm = sg.hm

    # Build dummy atoms
    dummy_atoms = []
    for w in selected:
        element = assignment[w.letter]
        coords = representative_coord(w.coordinates)
        dummy_atoms.append(DummyAtom(f"WYCK_{w.letter}", element, coords))

    # Export
    writers = {"cif": write_cif, "xyz": write_xyz}
    outputs = []
    for fmt in formats:
        if fmt in writers:
            path = f"{output_base}.{fmt}"
            writers[fmt](structure, dummy_atoms, path)
            outputs.append(path)

    return ", ".join(outputs)
