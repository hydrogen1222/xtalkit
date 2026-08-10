"""Core: mark Wyckoff positions in a crystal structure."""

import os
import math
import warnings

import gemmi

from xtalkit.spacegroup import wyckoff_positions
from xtalkit.matcher import match_atoms
from xtalkit.exporter import DummyAtom, write_cif, write_xyz
from xtalkit.builder import validate_cell
from xtalkit.utils import (
    assign_dummy_elements,
    read_cif_structure,
    representative_coord,
)


def mark(
    cif_path: str,
    sg_number: int,
    wyckoff_letters: list[str],
    mode: str,
    tolerance: float,
    element_map: dict[str, str] | None,
    formats: list[str],
    output_base: str,
    offset: float = 0.0,
) -> str:
    """Mark Wyckoff positions in a CIF file with dummy atoms.

    offset: optional visual displacement applied to dummy atoms.  A non-zero
            value deliberately moves them away from the exact Wyckoff orbit
            and is therefore unsuitable for crystallographic analysis.
    Returns comma-separated output file paths.
    """
    # Validate inputs
    if mode not in ("overlay", "replace"):
        raise ValueError(f"Mode must be 'overlay' or 'replace', got '{mode}'")
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be a finite non-negative number")
    if not math.isfinite(offset):
        raise ValueError("offset must be a finite number")
    if offset != 0.0:
        warnings.warn(
            "non-zero marker offset is visual-only and moves dummy atoms away "
            "from their exact Wyckoff positions",
            RuntimeWarning,
            stacklevel=2,
        )
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

    # Check file exists
    resolved = os.path.abspath(cif_path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"CIF file not found: {resolved}")

    # Read structure with proper CIF parsing + SG override
    structure = read_cif_structure(resolved, sg_number)
    cell_params = {
        "a": structure.cell.a,
        "b": structure.cell.b,
        "c": structure.cell.c,
        "alpha": structure.cell.alpha,
        "beta": structure.cell.beta,
        "gamma": structure.cell.gamma,
    }
    metric_issues = validate_cell(sg_number, cell_params)
    if metric_issues:
        raise ValueError(
            "CIF unit-cell metric is inconsistent with its space group: "
            + "; ".join(metric_issues)
        )

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

    # Overlay mode adds one generic representative for every requested Wyckoff
    # *type*. It is a visualization key, not a claim about a refined free
    # parameter. Replace mode is handled below and preserves the actual atom
    # coordinates and occupancies.
    dummy_atoms = []
    if mode == "overlay":
        for w in selected:
            element = assignment[w.letter]
            representative = representative_coord(w.coordinates)
            coords = tuple((value + offset) % 1.0 for value in representative)
            dummy_atoms.append(DummyAtom(f"WYCK_{w.letter}", element, coords))

    # Replace each matched atom in place. In particular, do not collapse all
    # atoms on one Wyckoff type into a single generic representative: that
    # would change the refined coordinates, occupancies, and disorder model.
    if mode == "replace":
        replaced_by_letter = {letter: 0 for letter in wyckoff_letters}
        reserved_labels = {
            atom.name
            for model in structure
            for chain in model
            for residue in chain
            for atom in residue
            if atom_mapping.get(atom.name) not in replaced_by_letter
        }
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        letter = atom_mapping.get(atom.name)
                        if letter not in replaced_by_letter:
                            continue
                        replaced_by_letter[letter] += 1
                        index = replaced_by_letter[letter]
                        atom.element = gemmi.Element(assignment[letter])
                        base_label = (
                            f"WYCK_{letter}" if index == 1
                            else f"WYCK_{letter}_{index}"
                        )
                        label = base_label
                        suffix = 2
                        while label in reserved_labels:
                            label = f"{base_label}_{suffix}"
                            suffix += 1
                        atom.name = label
                        reserved_labels.add(label)
                        if offset != 0.0:
                            frac = structure.cell.fractionalize(atom.pos)
                            shifted = gemmi.Fractional(
                                (frac.x + offset) % 1.0,
                                (frac.y + offset) % 1.0,
                                (frac.z + offset) % 1.0,
                            )
                            atom.pos = structure.cell.orthogonalize(shifted)
        unmatched = [
            letter for letter, count in replaced_by_letter.items() if count == 0
        ]
        if unmatched:
            warnings.warn(
                "replace mode found no matching atom for requested Wyckoff "
                f"position(s): {', '.join(unmatched)}",
                RuntimeWarning,
                stacklevel=2,
            )

    # Export
    writers = {"cif": write_cif, "xyz": write_xyz}
    outputs = []
    for fmt in formats:
        if fmt in writers:
            path = f"{output_base}.{fmt}"
            writers[fmt](structure, dummy_atoms, path)
            outputs.append(path)

    return ", ".join(outputs)
