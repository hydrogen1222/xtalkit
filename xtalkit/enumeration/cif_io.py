"""Small CIF reader/writer for SHRY preparation and verification."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, replace
from fractions import Fraction

import gemmi


@dataclass(frozen=True)
class AtomRow:
    label: str
    type_symbol: str
    x: str
    y: str
    z: str
    occupancy: str

    @property
    def frac_key(self) -> tuple[str, str, str]:
        return (_norm_coord(self.x), _norm_coord(self.y), _norm_coord(self.z))

    def with_occupancy(self, occ: Fraction) -> "AtomRow":
        return replace(self, occupancy=_format_fraction_decimal(occ))


@dataclass
class CifData:
    block_name: str
    cell: dict[str, str]
    spacegroup_name: str | None
    spacegroup_number: int | None
    symops: list[str]
    atoms: list[AtomRow]


def read_cif(path: str) -> CifData:
    """Read cell, symmetry, and atom-site rows from a standard CIF."""
    doc = gemmi.cif.read_file(path)
    block = doc.sole_block()
    cell = {
        "a": _required(block, "_cell_length_a"),
        "b": _required(block, "_cell_length_b"),
        "c": _required(block, "_cell_length_c"),
        "alpha": _required(block, "_cell_angle_alpha"),
        "beta": _required(block, "_cell_angle_beta"),
        "gamma": _required(block, "_cell_angle_gamma"),
    }
    sg_name = block.find_value("_symmetry_space_group_name_H-M")
    sg_number_raw = block.find_value("_symmetry_Int_Tables_number")
    sg_number = int(sg_number_raw) if sg_number_raw else None

    symops = list(block.find_values("_symmetry_equiv_pos_as_xyz"))
    if not symops:
        symops = list(block.find_values("_space_group_symop_operation_xyz"))
    symops = [_strip_quotes(op) for op in symops]
    if not symops and sg_name:
        try:
            symops = [op.triplet() for op in gemmi.SpaceGroup(_strip_quotes(sg_name)).operations()]
        except Exception:
            pass

    table = block.find("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    ])
    if len(table) == 0:
        raise ValueError("CIF has no _atom_site_ rows with fractional coordinates")
    atoms = [AtomRow(*[str(v) for v in row]) for row in table]
    return CifData(block.name, cell, sg_name, sg_number, symops, atoms)


def write_cif(data: CifData, path: str) -> None:
    """Write a SHRY-ready CIF preserving cell and symmetry operations."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    doc = gemmi.cif.Document()
    block = doc.add_new_block(data.block_name or "xtalkit_shry")
    block.set_pair("_cell_length_a", data.cell["a"])
    block.set_pair("_cell_length_b", data.cell["b"])
    block.set_pair("_cell_length_c", data.cell["c"])
    block.set_pair("_cell_angle_alpha", data.cell["alpha"])
    block.set_pair("_cell_angle_beta", data.cell["beta"])
    block.set_pair("_cell_angle_gamma", data.cell["gamma"])
    if data.spacegroup_name:
        block.set_pair("_symmetry_space_group_name_H-M", data.spacegroup_name)
    if data.spacegroup_number:
        block.set_pair("_symmetry_Int_Tables_number", str(data.spacegroup_number))
    if data.symops:
        loop = block.init_loop("_symmetry_equiv_pos_", ["site_id", "as_xyz"])
        for i, op in enumerate(data.symops, 1):
            loop.add_row([str(i), f"'{_strip_quotes(op)}'"])

    loop = block.init_loop("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    ])
    for atom in data.atoms:
        loop.add_row([
            atom.label,
            atom.type_symbol,
            atom.x,
            atom.y,
            atom.z,
            atom.occupancy,
        ])
    doc.write_file(path)


def copy_file(src: str, dst: str) -> None:
    parent = os.path.dirname(dst)
    if parent:
        os.makedirs(parent, exist_ok=True)
    shutil.copyfile(src, dst)


def _required(block, tag: str) -> str:
    value = block.find_value(tag)
    if value is None:
        raise ValueError(f"CIF missing required tag {tag}")
    return str(value)


def _strip_quotes(value: str) -> str:
    value = str(value).strip()
    if len(value) >= 2 and value[0] in "'\"" and value[-1] == value[0]:
        return value[1:-1]
    return value


def _norm_coord(value: str) -> str:
    return f"{float(value) % 1.0:.6f}"


def _format_fraction_decimal(value: Fraction) -> str:
    return f"{float(value):.8f}".rstrip("0").rstrip(".")
