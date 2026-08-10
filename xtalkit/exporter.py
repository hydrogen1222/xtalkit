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
    def w(value: float) -> float:
        return value - round(value)

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

    # ``write_cif`` is used only by marker/skeleton outputs. Dummy elements
    # are visualization proxies and must not be mistaken for a physical
    # composition in a downstream calculation.
    block.set_pair(
        "_audit_creation_method",
        "'xtalkit visualization-only Wyckoff markers; dummy species are not physical atoms'",
    )
    block.set_pair("_xtalkit_output_purpose", "visualization_only")

    # Cell parameters
    block.set_pair("_cell_length_a", str(cell.a))
    block.set_pair("_cell_length_b", str(cell.b))
    block.set_pair("_cell_length_c", str(cell.c))
    block.set_pair("_cell_angle_alpha", str(cell.alpha))
    block.set_pair("_cell_angle_beta", str(cell.beta))
    block.set_pair("_cell_angle_gamma", str(cell.gamma))

    # Space group (CIF requires quoting of values with spaces)
    sg_obj = None
    if structure.spacegroup_hm:
        block.set_pair(
            "_symmetry_space_group_name_H-M",
            f"'{structure.spacegroup_hm}'",
        )

    # Symmetry operations from the space group
        try:
            sg_obj = gemmi.SpaceGroup(structure.spacegroup_hm)
            block.set_pair("_symmetry_Int_Tables_number", str(sg_obj.number))
        except Exception:
            sg_obj = None
    if sg_obj is not None:
        sym_loop = block.init_loop("_symmetry_equiv_pos_", ["site_id", "as_xyz"])
        for i, op in enumerate(sg_obj.operations()):
            sym_loop.add_row([str(i + 1), f"'{op.triplet()}'"])

    # Standard small-molecule CIF atom loop.  In particular, occupancy must be
    # retained: marked disordered structures are still crystallographic data,
    # not merely pictures.
    loop = block.init_loop("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    ])

    # Original atom labels must remain unique. Duplicate labels make atom-site
    # identity ambiguous in CIF and also break replace-mode matching.
    used_labels: set[str] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name in used_labels:
                        raise ValueError(f"duplicate atom label {atom.name!r}")
                    used_labels.add(atom.name)
                    frac = cell.fractionalize(atom.pos)
                    loop.add_row([
                        atom.name,
                        atom.element.name,
                        f"{frac.x:.10f}",
                        f"{frac.y:.10f}",
                        f"{frac.z:.10f}",
                        # Gemmi stores occupancy as float32; seven significant
                        # digits preserve all meaningful input precision without
                        # exposing binary tails such as 0.6909999847.
                        f"{atom.occ:.7g}",
                    ])

    # Dummy atoms
    for da in dummy_atoms:
        label = da.label
        suffix = 2
        while label in used_labels:
            label = f"{da.label}_{suffix}"
            suffix += 1
        used_labels.add(label)
        loop.add_row([
            label,
            da.element,
            f"{da.fractional_coords[0]:.10f}",
            f"{da.fractional_coords[1]:.10f}",
            f"{da.fractional_coords[2]:.10f}",
            "1",
        ])

    doc.write_file(output_path)


def write_xyz(
    structure: gemmi.Structure,
    dummy_atoms: list[DummyAtom],
    output_path: str,
) -> None:
    """Write a symmetry-expanded visualization cell as an XYZ file.

    XYZ cannot store symmetry or occupancy, so asymmetric-unit atoms and
    Wyckoff representatives are explicitly expanded under the structure's
    space group. Partial-occupancy species are each shown on their full orbit;
    the result is therefore for visualization, not quantitative analysis.
    """
    _ensure_dir(output_path)

    all_atoms: list[tuple[str, float, float, float]] = []
    cell = structure.cell
    sg = gemmi.SpaceGroup(structure.spacegroup_hm) \
        if structure.spacegroup_hm else gemmi.SpaceGroup(1)

    # Expand original asymmetric-unit atoms.
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    for image in _expanded_orbit(frac, sg):
                        orth = cell.orthogonalize(image)
                        all_atoms.append((
                            atom.element.name, orth.x, orth.y, orth.z,
                        ))

    # Expand each dummy Wyckoff representative.
    for da in dummy_atoms:
        frac = gemmi.Fractional(*da.fractional_coords)
        for image in _expanded_orbit(frac, sg):
            orth = cell.orthogonalize(image)
            all_atoms.append((da.element, orth.x, orth.y, orth.z))

    lines = [
        f"{len(all_atoms)}",
        "xtalkit symmetry-expanded visualization; occupancy is not encoded",
    ]
    for el, x, y, z in all_atoms:
        lines.append(f"{el:<4} {x:12.6f} {y:12.6f} {z:12.6f}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def _expanded_orbit(
    representative: gemmi.Fractional,
    spacegroup: gemmi.SpaceGroup,
) -> list[gemmi.Fractional]:
    """Expand and periodically deduplicate one fractional representative."""
    images: list[gemmi.Fractional] = []
    for op in spacegroup.operations():
        transformed = op.apply_to_xyz([
            representative.x, representative.y, representative.z,
        ])
        image = gemmi.Fractional(*transformed)
        image.wrap_to_unit()
        if not any(_frac_dist(image, old) < 1e-7 for old in images):
            images.append(image)
    return images


def _chemical_formula_sum(
    structure: gemmi.Structure,
    spacegroup: gemmi.SpaceGroup,
) -> str:
    """Return the symmetry-expanded unit-cell composition for a CIF tag."""
    amounts: dict[str, float] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = structure.cell.fractionalize(atom.pos)
                    multiplicity = len(_expanded_orbit(frac, spacegroup))
                    amounts[atom.element.name] = (
                        amounts.get(atom.element.name, 0.0)
                        + float(atom.occ) * multiplicity
                    )
    if not amounts:
        raise ValueError("cannot write a CIF for a structure without atoms")
    parts = []
    for symbol, amount in amounts.items():
        normalized = round(amount, 7)
        if normalized <= 0:
            raise ValueError(
                f"non-positive unit-cell amount for {symbol}: {normalized}"
            )
        parts.append(f"{symbol}{normalized:.7g}")
    return " ".join(parts)


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

    formula_group = sg or gemmi.SpaceGroup(1)
    block.set_pair(
        "_chemical_formula_sum",
        f"'{_chemical_formula_sum(structure, formula_group)}'",
    )

    loop = block.init_loop("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    ])
    used_labels: set[str] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name in used_labels:
                        raise ValueError(f"duplicate atom label {atom.name!r}")
                    used_labels.add(atom.name)
                    frac = cell.fractionalize(atom.pos)
                    loop.add_row([
                        atom.name,
                        atom.element.name,
                        f"{frac.x:.10f}",
                        f"{frac.y:.10f}",
                        f"{frac.z:.10f}",
                        f"{atom.occ:.7g}",
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
                    for f in _expanded_orbit(frac, sg):
                        orth = cell.orthogonalize(f)
                        all_atoms.append(
                            (atom.element.name, orth.x, orth.y, orth.z))

    lines = [
        str(len(all_atoms)),
        "xtalkit symmetry-expanded visualization; occupancy is not encoded",
    ]
    for el, x, y, z in all_atoms:
        lines.append(f"{el:<4} {x:12.6f} {y:12.6f} {z:12.6f}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))
