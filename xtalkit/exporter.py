"""Output writers for .cif and .xyz formats."""

import math
import os
from collections import namedtuple

import gemmi

DummyAtom = namedtuple("DummyAtom", ["label", "element", "fractional_coords"])


def _ensure_dir(path: str) -> None:
    """Ensure the parent directory exists."""
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def _frac_dist(a: gemmi.Fractional, b: gemmi.Fractional) -> float:
    """Minimum-image distance in fractional space (for orbit dedup)."""
    w = lambda x: x - round(x)
    return math.sqrt(w(a.x - b.x) ** 2 + w(a.y - b.y) ** 2 + w(a.z - b.z) ** 2)


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

    # Symmetry operations from the space group
    try:
        sg_obj = gemmi.SpaceGroup(structure.spacegroup_hm)
    except Exception:
        sg_obj = None
    if sg_obj is not None:
        sym_loop = block.init_loop("_symmetry_equiv_pos_", ["site_id", "as_xyz"])
        for i, op in enumerate(sg_obj.operations()):
            sym_loop.add_row([str(i + 1), f"'{op.triplet()}'"])

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


def write_structure_cif(
    structure: gemmi.Structure,
    output_path: str,
) -> None:
    """Write an asymmetric-unit CIF (standard crystallographic tags).

    Emits ``_cell_*``, ``_symmetry_space_group_name_H-M``, the
    ``_symmetry_equiv_pos_as_xyz`` loop, and an ``_atom_site_`` loop with
    fractional coordinates and occupancy. Only the asymmetric-unit
    representatives are written; the reader expands them via the symmetry
    operations. This is the format VESTA, GSAS-II, and pymatgen expect for a
    refinement CIF.
    """
    _ensure_dir(output_path)

    cell = structure.cell
    doc = gemmi.cif.Document()
    block = doc.add_new_block(structure.name or "xtalkit_built")

    block.set_pair("_cell_length_a", f"{cell.a:.6f}")
    block.set_pair("_cell_length_b", f"{cell.b:.6f}")
    block.set_pair("_cell_length_c", f"{cell.c:.6f}")
    block.set_pair("_cell_angle_alpha", f"{cell.alpha:.6f}")
    block.set_pair("_cell_angle_beta", f"{cell.beta:.6f}")
    block.set_pair("_cell_angle_gamma", f"{cell.gamma:.6f}")

    sg = None
    if structure.spacegroup_hm:
        block.set_pair(
            "_symmetry_space_group_name_H-M",
            f"'{structure.spacegroup_hm}'",
        )
        try:
            sg = gemmi.SpaceGroup(structure.spacegroup_hm)
            block.set_pair("_symmetry_Int_Tables_number", str(sg.number))
        except Exception:
            sg = None

    if sg is not None:
        sym_loop = block.init_loop(
            "_symmetry_equiv_pos_", ["site_id", "as_xyz"])
        for i, op in enumerate(sg.operations()):
            sym_loop.add_row([str(i + 1), f"'{op.triplet()}'"])

    loop = block.init_loop("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    ])
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    loop.add_row([
                        atom.name,
                        atom.element.name,
                        f"{frac.x:.6f}",
                        f"{frac.y:.6f}",
                        f"{frac.z:.6f}",
                        f"{atom.occ:.4f}",
                    ])

    doc.write_file(output_path)


def write_structure_xyz(
    structure: gemmi.Structure,
    sg_number: int,
    output_path: str,
) -> None:
    """Write a symmetry-expanded XYZ (full unit cell, Cartesian).

    Each asymmetric-unit representative is expanded under the space group's
    operations (deduplicated) so the XYZ shows every atom in the cell. XYZ
    cannot encode occupancy, so partial-occupancy sites emit each species'
    full orbit — the file is meant for quick visualization, not analysis.
    """
    _ensure_dir(output_path)

    cell = structure.cell
    sg = gemmi.SpaceGroup(sg_number)

    all_atoms: list[tuple[str, float, float, float]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    seen: list[gemmi.Fractional] = []
                    for op in sg.operations():
                        c = op.apply_to_xyz([frac.x, frac.y, frac.z])
                        f = gemmi.Fractional(*c)
                        f.wrap_to_unit()
                        if not any(_frac_dist(f, s) < 1e-3 for s in seen):
                            seen.append(f)
                    for f in seen:
                        orth = cell.orthogonalize(f)
                        all_atoms.append(
                            (atom.element.name, orth.x, orth.y, orth.z))

    lines = [str(len(all_atoms)), "xtalkit built"]
    for el, x, y, z in all_atoms:
        lines.append(f"{el:<4} {x:12.6f} {y:12.6f} {z:12.6f}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))
