"""Generate pure Wyckoff skeleton structures (no real atoms)."""

import os

import gemmi
from xtalkit.spacegroup import wyckoff_positions, default_cell_params
from xtalkit.exporter import DummyAtom, write_cif, write_vesta, write_xyz


_DUMMY_ELEMENTS = ["Xe", "Kr", "Rn", "Ar", "Ne", "He"]
_VARIABLE_DEFAULT = 0.125


def _assign_dummy_elements(
    wyckoff_letters: list[str],
    element_map: dict[str, str] | None,
) -> dict[str, str]:
    """Assign dummy elements to Wyckoff letters."""
    if element_map:
        return dict(element_map)
    assignment = {}
    for i, letter in enumerate(sorted(
        wyckoff_letters,
        key=lambda w: (int("".join(c for c in w if c.isdigit()) or 0), w),
    )):
        assignment[letter] = _DUMMY_ELEMENTS[i % len(_DUMMY_ELEMENTS)]
    return assignment


def _parse_coord(s: str) -> float:
    """Parse a coordinate expression like '0', '1/4', '0.25', 'x' to float."""
    s = s.strip()
    if s in ("x", "y", "z"):
        return _VARIABLE_DEFAULT
    if "/" in s:
        num, den = s.split("/")
        return float(num) / float(den)
    return float(s)


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
    assignment = _assign_dummy_elements(wyckoff_letters, element_map)

    # Cell parameters
    if cell_params is None:
        params = default_cell_params(sg_number)
    else:
        params = cell_params

    # Build structure
    structure = gemmi.Structure()
    structure.cell = gemmi.UnitCell(
        params["a"], params["b"], params["c"],
        params["alpha"], params["beta"], params["gamma"],
    )
    sg = gemmi.SpaceGroup(sg_number)
    structure.spacegroup_hm = sg.hm

    # Build dummy atoms
    dummy_atoms = []
    for w in selected:
        element = assignment[w.letter]
        parts = w.coordinates.split(",")
        coords = tuple(_parse_coord(p.strip()) for p in parts)
        dummy_atoms.append(DummyAtom(f"WYCK_{w.letter}", element, coords))

    # Place them in structure (as a single model/chain/residue)
    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "WYC"
    residue.seqid = gemmi.SeqId("1")

    for da in dummy_atoms:
        atom = gemmi.Atom()
        atom.name = da.label
        atom.element = gemmi.Element(da.element)
        # gemmi 0.7.5: use Cartesian coordinates via orthogonalize
        frac = gemmi.Fractional(*da.fractional_coords)
        orth = structure.cell.orthogonalize(frac)
        atom.pos = (orth.x, orth.y, orth.z)
        residue.add_atom(atom)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)

    # Export
    writers = {"cif": write_cif, "vesta": write_vesta, "xyz": write_xyz}
    outputs = []
    for fmt in formats:
        if fmt in writers:
            path = f"{output_base}.{fmt}"
            writers[fmt](structure, dummy_atoms, path)
            outputs.append(path)

    return ", ".join(outputs)
