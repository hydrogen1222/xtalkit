"""Core: mark Wyckoff positions in a crystal structure."""

import os

import gemmi
from xtalkit.spacegroup import wyckoff_positions
from xtalkit.matcher import match_atoms
from xtalkit.exporter import DummyAtom, write_cif, write_vesta, write_xyz
from xtalkit.utils import assign_dummy_elements, parse_coord


def mark(
    cif_path: str,
    sg_number: int,
    wyckoff_letters: list[str],
    mode: str,
    tolerance: float,
    element_map: dict[str, str] | None,
    formats: list[str],
    output_base: str,
) -> str:
    """Mark Wyckoff positions in a CIF file with dummy atoms.

    Returns comma-separated output file paths.
    """
    # Validate inputs
    if mode not in ("overlay", "replace"):
        raise ValueError(f"Mode must be 'overlay' or 'replace', got '{mode}'")

    # Check file exists
    resolved = os.path.abspath(cif_path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"CIF file not found: {resolved}")

    # Read structure
    structure = gemmi.read_structure(resolved)

    # Get Wyckoff positions for this SG
    all_wyckoffs = wyckoff_positions(sg_number)
    valid_letters = {w.letter for w in all_wyckoffs}

    # Validate requested letters
    for letter in wyckoff_letters:
        if letter not in valid_letters:
            raise ValueError(
                f"Invalid Wyckoff letter '{letter}' for SG #{sg_number}. "
                f"Valid letters: {', '.join(sorted(valid_letters))}"
            )

    # Filter to requested Wyckoff positions
    selected = [w for w in all_wyckoffs if w.letter in wyckoff_letters]

    # Assign dummy elements
    assignment = assign_dummy_elements(wyckoff_letters, element_map)

    # Match atoms if needed for replace mode
    atom_mapping = {}
    if mode == "replace":
        atom_mapping = match_atoms(structure, all_wyckoffs, tolerance)

    # Build dummy atom list
    dummy_atoms = []
    for w in selected:
        element = assignment[w.letter]
        # Parse coordinates string to tuple
        parts = w.coordinates.split(",")
        coords = tuple(parse_coord(p.strip()) for p in parts)
        dummy_atoms.append(DummyAtom(f"WYCK_{w.letter}", element, coords))

    # Handle replace mode: remove original atoms at matched positions,
    # replace with dummy atoms
    if mode == "replace":
        labels_to_remove = {
            label for label, letter in atom_mapping.items()
            if letter in wyckoff_letters
        }
        for model in structure:
            for chain in model:
                for residue in chain:
                    to_remove = [a for a in residue if a.name in labels_to_remove]
                    for atom in to_remove:
                        residue.remove_atom(atom.name, ".", atom.element)  # type: ignore

    # Export
    writers = {"cif": write_cif, "vesta": write_vesta, "xyz": write_xyz}
    outputs = []
    for fmt in formats:
        if fmt in writers:
            path = f"{output_base}.{fmt}"
            writers[fmt](structure, dummy_atoms, path)
            outputs.append(path)

    return ", ".join(outputs)
