"""Output writers for .cif, .vesta, and .xyz formats."""

import os
from collections import namedtuple

import gemmi

DummyAtom = namedtuple("DummyAtom", ["label", "element", "fractional_coords"])


def _ensure_dir(path: str) -> None:
    """Ensure the parent directory exists."""
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def write_cif(
    structure: gemmi.Structure,
    dummy_atoms: list[DummyAtom],
    output_path: str,
) -> None:
    """Write structure with dummy atoms as a CIF file.

    Both original atoms and dummy atoms are included with fractional coordinates.
    """
    _ensure_dir(output_path)

    cell = structure.cell
    doc = gemmi.cif.Document()
    block_name = structure.name if structure.name else "xtalkit"
    block = doc.add_new_block(block_name)

    # Cell parameters
    block.set_pair("_cell_length_a", str(cell.a))
    block.set_pair("_cell_length_b", str(cell.b))
    block.set_pair("_cell_length_c", str(cell.c))
    block.set_pair("_cell_angle_alpha", str(cell.alpha))
    block.set_pair("_cell_angle_beta", str(cell.beta))
    block.set_pair("_cell_angle_gamma", str(cell.gamma))

    # Space group (CIF requires quoting of values with spaces)
    if structure.spacegroup_hm:
        block.set_pair(
            "_symmetry_space_group_name_H-M",
            f"'{structure.spacegroup_hm}'",
        )

    # Atom site loop (mmCIF format with fractional coords)
    loop = block.init_loop("_atom_site.", [
        "group_PDB", "id", "type_symbol", "label_atom_id",
        "Cartn_x", "Cartn_y", "Cartn_z",
        "fract_x", "fract_y", "fract_z",
    ])

    # Original atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    loop.add_row([
                        "ATOM",
                        atom.name,
                        atom.element.name,
                        atom.name,
                        f"{atom.pos.x:.6f}",
                        f"{atom.pos.y:.6f}",
                        f"{atom.pos.z:.6f}",
                        f"{frac.x:.6f}",
                        f"{frac.y:.6f}",
                        f"{frac.z:.6f}",
                    ])

    # Dummy atoms
    for da in dummy_atoms:
        frac = gemmi.Fractional(*da.fractional_coords)
        orth = cell.orthogonalize(frac)
        loop.add_row([
            "ATOM",
            da.label,
            da.element,
            da.label,
            f"{orth.x:.6f}",
            f"{orth.y:.6f}",
            f"{orth.z:.6f}",
            f"{da.fractional_coords[0]:.6f}",
            f"{da.fractional_coords[1]:.6f}",
            f"{da.fractional_coords[2]:.6f}",
        ])

    doc.write_file(output_path)


def write_vesta(
    structure: gemmi.Structure,
    dummy_atoms: list[DummyAtom],
    output_path: str,
) -> None:
    """Write structure with dummy atoms as a VESTA (.vesta) file.

    VESTA format is XML-based. We generate a minimal structure document.
    """
    _ensure_dir(output_path)

    cell = structure.cell
    all_atoms: list[tuple[str, float, float, float]] = []

    # Collect original atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    all_atoms.append((atom.element.name, frac.x, frac.y, frac.z))

    # Add dummy atoms
    for da in dummy_atoms:
        all_atoms.append((da.element, *da.fractional_coords))

    # Generate VESTA XML
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<VESTA xmlns="http://www.jp-minerals.org/vesta/">',
        "  <data>",
        "    <structure>",
        f'      <unit_cell a="{cell.a}" b="{cell.b}" c="{cell.c}"',
        f'                 alpha="{cell.alpha}" beta="{cell.beta}" gamma="{cell.gamma}"/>',
        "      <site_list>",
    ]

    for i, (el, x, y, z) in enumerate(all_atoms):
        lines.extend([
            "        <site>",
            f"          <element>{el}</element>",
            f"          <x>{x:.8f}</x>",
            f"          <y>{y:.8f}</y>",
            f"          <z>{z:.8f}</z>",
            f"          <label>{el}{i+1}</label>",
            "        </site>",
        ])

    lines.extend([
        "      </site_list>",
        "    </structure>",
        "  </data>",
        "</VESTA>",
    ])

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def write_xyz(
    structure: gemmi.Structure,
    dummy_atoms: list[DummyAtom],
    output_path: str,
) -> None:
    """Write structure with dummy atoms as an XYZ file.

    Coordinates are converted to Cartesian using the structure's cell.
    """
    _ensure_dir(output_path)

    all_atoms: list[tuple[str, float, float, float]] = []
    cell = structure.cell

    # Collect original atoms (atom.pos is already Cartesian)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    all_atoms.append((
                        atom.element.name,
                        atom.pos.x,
                        atom.pos.y,
                        atom.pos.z,
                    ))

    # Add dummy atoms (convert fractional to Cartesian)
    for da in dummy_atoms:
        frac = gemmi.Fractional(*da.fractional_coords)
        orth = cell.orthogonalize(frac)
        all_atoms.append((da.element, orth.x, orth.y, orth.z))

    lines = [f"{len(all_atoms)}", "xtalkit generated"]
    for el, x, y, z in all_atoms:
        lines.append(f"{el:<4} {x:12.6f} {y:12.6f} {z:12.6f}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))
