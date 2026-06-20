"""Shared utility functions for xtalkit."""

import re

import gemmi

_DUMMY_ELEMENTS = ["Xe", "Kr", "Rn", "Ar", "Ne", "He"]


def read_cif_structure(path: str, sg_number: int) -> gemmi.Structure:
    """Read a CIF file and build a gemmi.Structure with the given space group.

    Uses gemmi.cif.read_file for robust parsing. Handles both standard CIF
    (_cell_length_a, _atom_site_fract_x) and mmCIF (_cell.length_a,
    _atom_site.Cartn_x) tag conventions. Falls back to gemmi.read_structure
    if neither format is detected.
    """
    doc = gemmi.cif.read_file(path)
    block = doc.sole_block()

    # Try standard CIF cell tags first, then mmCIF
    def _cell_val(*names: str) -> float | None:
        for name in names:
            v = block.find_value(name)
            if v is not None:
                return float(v)
        return None

    a = _cell_val("_cell_length_a", "_cell.length_a")
    b = _cell_val("_cell_length_b", "_cell.length_b")
    c = _cell_val("_cell_length_c", "_cell.length_c")
    alpha = _cell_val("_cell_angle_alpha", "_cell.angle_alpha")
    beta = _cell_val("_cell_angle_beta", "_cell.angle_beta")
    gamma = _cell_val("_cell_angle_gamma", "_cell.angle_gamma")

    if None in (a, b, c, alpha, beta, gamma):
        raise ValueError(f"Could not read cell parameters from CIF: {path}")

    cell = gemmi.UnitCell(a, b, c, alpha, beta, gamma)

    # Try to read atoms from standard CIF loop first, then mmCIF
    labels = list(block.find_values("_atom_site_label"))
    if not labels:
        labels = list(block.find_values("_atom_site.label_atom_id"))

    type_col = list(block.find_values("_atom_site_type_symbol"))
    if not type_col:
        type_col = list(block.find_values("_atom_site.type_symbol"))

    fx = list(block.find_values("_atom_site_fract_x"))
    fy = list(block.find_values("_atom_site_fract_y"))
    fz = list(block.find_values("_atom_site_fract_z"))

    # mmCIF uses Cartesian; convert to fractional
    cart_x = list(block.find_values("_atom_site.Cartn_x"))
    cart_y = list(block.find_values("_atom_site.Cartn_y"))
    cart_z = list(block.find_values("_atom_site.Cartn_z"))

    use_cartesian = False
    if not fx:
        fx, fy, fz = cart_x, cart_y, cart_z
        use_cartesian = True

    n_atoms = len(labels)
    if n_atoms == 0:
        raise ValueError("No atom sites found in CIF file")

    # Build structure
    structure = gemmi.Structure()
    structure.cell = cell
    sg = gemmi.SpaceGroup(sg_number)
    structure.spacegroup_hm = sg.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    for i in range(n_atoms):
        raw_element = type_col[i] if i < len(type_col) else "?"
        element = re.sub(r"[0-9+\-]", "", raw_element) if raw_element else "?"
        if not element or element == "?":
            element = "Xe"

        atom = gemmi.Atom()
        atom.name = labels[i] if i < len(labels) else f"X{i}"
        atom.element = gemmi.Element(element)

        if use_cartesian and i < len(fx):
            atom.pos = gemmi.Position(float(fx[i]), float(fy[i]), float(fz[i]))
        elif i < len(fx):
            frac = gemmi.Fractional(float(fx[i]), float(fy[i]), float(fz[i]))
            atom.pos = cell.orthogonalize(frac)
        else:
            atom.pos = gemmi.Position(0.0, 0.0, 0.0)

        residue.add_atom(atom)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)

    return structure


def assign_dummy_elements(
    wyckoff_letters: list[str],
    element_map: dict[str, str] | None,
) -> dict[str, str]:
    """Assign dummy elements to Wyckoff letters.

    Priority: Xe -> Kr -> Rn -> Ar -> Ne -> He (cycling if needed).

    When element_map is provided, every requested letter must have an entry.
    """
    if element_map is not None:
        for letter in wyckoff_letters:
            if letter not in element_map:
                raise ValueError(f"No element assigned for Wyckoff letter {letter}")
        return dict(element_map)

    assignment = {}
    for i, letter in enumerate(sorted(
        wyckoff_letters,
        key=lambda w: (int("".join(c for c in w if c.isdigit()) or 0), w),
    )):
        assignment[letter] = _DUMMY_ELEMENTS[i % len(_DUMMY_ELEMENTS)]
    return assignment


def parse_coord(s: str, variable_default: float = 0.3) -> float:
    """Parse a coordinate expression like '0', '1/4', '0.25' to float.

    Variable expressions ('x', 'y', 'z') are resolved to a representative
    value (default 0.3) for dummy atom placement.
    """
    s = s.strip()
    if s in ("x", "y", "z"):
        return variable_default
    if "/" in s:
        num, den = s.split("/")
        return float(num) / float(den)
    return float(s)
