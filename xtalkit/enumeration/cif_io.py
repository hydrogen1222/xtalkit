"""Small CIF reader/writer for SHRY preparation and verification."""

from __future__ import annotations

import os
import math
import re
import shutil
from dataclasses import dataclass, replace
from fractions import Fraction

import gemmi

from xtalkit.utils import unit_cell_from_parameters


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
        "a": _strip_esd(_required(block, "_cell_length_a")),
        "b": _strip_esd(_required(block, "_cell_length_b")),
        "c": _strip_esd(_required(block, "_cell_length_c")),
        "alpha": _strip_esd(_required(block, "_cell_angle_alpha")),
        "beta": _strip_esd(_required(block, "_cell_angle_beta")),
        "gamma": _strip_esd(_required(block, "_cell_angle_gamma")),
    }
    try:
        unit_cell_from_parameters(cell)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid CIF unit cell in {path}: {exc}") from exc
    sg_name = (
        block.find_value("_symmetry_space_group_name_H-M")
        or block.find_value("_space_group_name_H-M_alt")
    )
    sg_number_raw = (
        block.find_value("_symmetry_Int_Tables_number")
        or block.find_value("_space_group_IT_number")
    )
    try:
        sg_number = int(sg_number_raw) if sg_number_raw else None
    except ValueError as exc:
        raise ValueError(
            f"invalid CIF International Tables number {sg_number_raw!r}"
        ) from exc
    if sg_number is not None and not 1 <= sg_number <= 230:
        raise ValueError(f"CIF space-group number must be 1-230, got {sg_number}")

    declared_group = None
    if sg_name:
        try:
            declared_group = gemmi.SpaceGroup(_strip_quotes(sg_name))
        except Exception as exc:
            raise ValueError(f"invalid CIF space-group name {sg_name!r}") from exc
        if sg_number is not None and declared_group.number != sg_number:
            raise ValueError(
                "CIF space-group name and International Tables number disagree: "
                f"{_strip_quotes(sg_name)!r} is #{declared_group.number}, "
                f"not #{sg_number}"
            )
    elif sg_number is not None:
        declared_group = gemmi.find_spacegroup_by_number(sg_number)
    if sg_number is None and declared_group is not None:
        sg_number = declared_group.number

    symops = list(block.find_values("_symmetry_equiv_pos_as_xyz"))
    if not symops:
        symops = list(block.find_values("_space_group_symop_operation_xyz"))
    symops = [_strip_quotes(op) for op in symops]
    if not symops and declared_group is not None:
        symops = [op.triplet() for op in declared_group.operations()]

    labels = list(block.find_values("_atom_site_label"))
    types = list(block.find_values("_atom_site_type_symbol"))
    xs = list(block.find_values("_atom_site_fract_x"))
    ys = list(block.find_values("_atom_site_fract_y"))
    zs = list(block.find_values("_atom_site_fract_z"))
    occupancies = list(block.find_values("_atom_site_occupancy"))
    if not labels:
        raise ValueError("CIF has no _atom_site_ rows with fractional coordinates")
    n_rows = len(labels)
    if len(set(labels)) != n_rows:
        raise ValueError("CIF atom-site labels must be unique")
    if not (len(xs) == len(ys) == len(zs) == n_rows):
        raise ValueError(
            "CIF atom-site coordinate columns have inconsistent lengths: "
            f"labels/x/y/z={n_rows}/{len(xs)}/{len(ys)}/{len(zs)}"
        )
    if types and len(types) != n_rows:
        raise ValueError(
            f"CIF atom-site type_symbol has {len(types)} rows, expected {n_rows}"
        )
    if occupancies and len(occupancies) != n_rows:
        raise ValueError(
            f"CIF atom-site occupancy has {len(occupancies)} rows, expected {n_rows}"
        )
    atoms = []
    for index, label in enumerate(labels):
        type_symbol = types[index] if types else _element_from_label(label)
        occ = occupancies[index] if occupancies else "1"
        coordinate_values = [
            _strip_esd(xs[index]), _strip_esd(ys[index]), _strip_esd(zs[index]),
        ]
        try:
            numeric_coords = [float(value) for value in coordinate_values]
        except ValueError as exc:
            raise ValueError(
                f"invalid fractional coordinate for atom {label!r}"
            ) from exc
        if not all(math.isfinite(value) for value in numeric_coords):
            raise ValueError(f"non-finite fractional coordinate for atom {label!r}")
        try:
            Fraction(_strip_esd(occ))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"invalid occupancy {occ!r} for atom {label!r}") from exc
        atoms.append(AtomRow(
            str(label), str(type_symbol),
            *coordinate_values, _strip_esd(occ),
        ))
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
            _occupancy_decimal(atom.occupancy),
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


def _element_from_label(label: str) -> str:
    match = re.match(r"[A-Z][a-z]?", str(label).strip())
    if not match:
        raise ValueError(
            "CIF is missing _atom_site_type_symbol and element cannot be "
            f"inferred from label {label!r}"
        )
    return match.group(0)


def _norm_coord(value: str) -> str:
    rounded = round(float(value) % 1.0, 6) % 1.0
    return f"{rounded:.6f}"


_ESD_RE = re.compile(r"(?<=\d)\(\d+\)(?=(?:[eEdD][+-]?\d+)?$)")


def _strip_esd(value: str) -> str:
    """Strip a trailing CIF standard uncertainty, e.g. ``0.3148(19)`` -> ``0.3148``.

    Real-world CIFs from refinements (ICSD, VESTA exports) carry esds on cell
    lengths and fractional coordinates; ``float()`` chokes on them. Labels and
    type symbols are passed through unchanged.
    """
    s = str(value).strip()
    return _ESD_RE.sub("", s).replace("D", "E").replace("d", "e")


def _format_fraction_decimal(value: Fraction) -> str:
    return f"{float(value):.16g}"


def _occupancy_decimal(value: str) -> str:
    """Normalize an occupancy string to a CIF-safe decimal.

    SHRY parses CIFs with pymatgen, whose ``float()``-based occupancy reader
    rejects fraction strings like ``1/2``. Inputs may carry fractions (from
    :func:`fraction_string`), decimals (``0.5``), or ``1``; emit a plain
    decimal so downstream parsers stay happy.
    """
    return _format_fraction_decimal(Fraction(str(value).strip()))
