"""Build SHRY-ready CIFs with strict occupancy and symmetry checks."""

from __future__ import annotations

import os
from fractions import Fraction

from xtalkit.enumeration.cif_io import AtomRow, read_cif, write_cif
from xtalkit.enumeration.manifest import base_manifest, sha256_file, write_json
from xtalkit.enumeration.occupancy import (
    check_integerizable,
    fraction_string,
    parse_fraction,
    parse_set_occupancy,
)
from xtalkit.enumeration.symmetry_audit import audit_symmetry_grouping


def prepare_shry_input(
    input_cif: str,
    output_cif: str,
    vacancy_symbol: str = "X",
    occupancy_overrides: str | None = None,
    parent_spacegroup: int | None = None,
    target_formula: str | None = None,
    strict: bool = True,
    symmetrize: bool = False,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
    scaling_matrix: list[list[int]] | None = None,
) -> dict:
    """Prepare a SHRY-ready CIF and sidecar manifest."""
    if not vacancy_symbol:
        raise ValueError("vacancy_symbol must be non-empty")
    if not os.path.exists(input_cif):
        raise FileNotFoundError(input_cif)

    matrix = scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    det = _det3(matrix)
    overrides = parse_set_occupancy(occupancy_overrides)
    data = read_cif(input_cif)
    data.atoms = _apply_overrides(data.atoms, overrides)
    data.atoms = _fill_and_validate_vacancies(data.atoms, vacancy_symbol)

    audit = audit_symmetry_grouping(
        data,
        expected_spacegroup=parent_spacegroup,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        symmetrize=symmetrize,
    )
    _check_integerizable_groups(data.atoms, audit["site_orbits"], det)

    write_cif(data, output_cif)
    manifest_path = f"{output_cif}.manifest.json"
    orbit_path = f"{output_cif}.orbit_grouping.json"
    manifest = base_manifest(
        mode="strict_exhaustive" if strict else "shry_prepare",
        input_cif=os.path.abspath(output_cif),
        original_cif=os.path.abspath(input_cif),
        input_sha256=sha256_file(input_cif),
        prepared_sha256=sha256_file(output_cif),
        parent_spacegroup_declared=parent_spacegroup,
        parent_spacegroup_detected=audit["parent_spacegroup_detected"],
        hall_symbol=audit["hall_symbol"],
        site_orbits=audit["site_orbits"],
        occupancy_idealization=_override_manifest(overrides),
        supercell_matrix=matrix,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        symmetrize_flag=bool(symmetrize),
        vacancy_symbol=vacancy_symbol,
        target_formula_after_removing_vacancy=target_formula,
    )
    write_json(manifest_path, manifest)
    write_json(orbit_path, audit)
    return {
        "output_cif": output_cif,
        "manifest": manifest_path,
        "orbit_grouping": orbit_path,
        "site_orbits": audit["site_orbits"],
    }


def _apply_overrides(
    atoms: list[AtomRow],
    overrides: dict[str, dict[str, Fraction]],
) -> list[AtomRow]:
    if not overrides:
        return atoms
    consumed: set[str] = set()
    out: list[AtomRow] = []
    for atom in atoms:
        override = overrides.get(atom.label)
        if override is None:
            out.append(atom)
            continue
        consumed.add(atom.label)
        if "__self__" in override:
            out.append(atom.with_occupancy(override["__self__"]))
        else:
            for species, occ in override.items():
                out.append(AtomRow(
                    label=f"{atom.label}_{species}",
                    type_symbol=species,
                    x=atom.x,
                    y=atom.y,
                    z=atom.z,
                    occupancy=fraction_string(occ),
                ))
    unknown = set(overrides) - consumed
    if unknown:
        raise ValueError(f"--set-occupancy label(s) not found in CIF: {sorted(unknown)}")
    return out


def _fill_and_validate_vacancies(
    atoms: list[AtomRow],
    vacancy_symbol: str,
) -> list[AtomRow]:
    groups: dict[tuple[str, str, str], list[AtomRow]] = {}
    for atom in atoms:
        groups.setdefault(atom.frac_key, []).append(atom)

    out: list[AtomRow] = []
    for group in groups.values():
        total = sum(parse_fraction(atom.occupancy) for atom in group)
        if total > 1:
            labels = ", ".join(atom.label for atom in group)
            raise ValueError(f"occupancy sum exceeds 1 for site {labels}: {total}")
        out.extend(group)
        if total < 1:
            base = group[0]
            missing = 1 - total
            out.append(AtomRow(
                label=f"{base.label}_{vacancy_symbol}",
                type_symbol=vacancy_symbol,
                x=base.x,
                y=base.y,
                z=base.z,
                occupancy=fraction_string(missing),
            ))
    return out


def _check_integerizable_groups(
    atoms: list[AtomRow],
    site_orbits: list[dict],
    det: int,
) -> None:
    atoms_by_coord = {atom.frac_key: [] for atom in atoms}
    for atom in atoms:
        atoms_by_coord[atom.frac_key].append(atom)

    for record in site_orbits:
        rep = tuple(record["representative"])
        group_atoms = atoms_by_coord.get(rep, [])
        if not group_atoms:
            continue
        occupancies = {
            atom.type_symbol: parse_fraction(atom.occupancy)
            for atom in group_atoms
        }
        n_sites = int(record["multiplicity"]) * det
        check_integerizable(n_sites, occupancies, ",".join(record["labels"]))
        total = sum(occupancies.values())
        if total != 1:
            raise ValueError(
                f"occupancy sum for orbit {record['orbit_id']} is "
                f"{fraction_string(total)}, expected 1"
            )


def _override_manifest(overrides: dict[str, dict[str, Fraction]]) -> dict:
    return {
        key: {sp: fraction_string(occ) for sp, occ in value.items()}
        for key, value in overrides.items()
    }


def _det3(matrix: list[list[int]]) -> int:
    a = matrix
    det = (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )
    if det <= 0:
        raise ValueError(f"scaling matrix determinant must be positive, got {det}")
    return det
