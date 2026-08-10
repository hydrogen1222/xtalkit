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
from xtalkit.enumeration.formula import (
    element_symbol,
    parse_formula,
    real_element_symbol,
    validate_vacancy_symbol,
)
from xtalkit.enumeration.symmetry_audit import audit_symmetry_grouping
from xtalkit.enumeration.shry_backend import (
    scaling_determinant,
    validate_scaling_matrix,
    validate_tolerances,
)


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
    vacancy_symbol = validate_vacancy_symbol(vacancy_symbol)
    if not os.path.exists(input_cif):
        raise FileNotFoundError(input_cif)
    validate_tolerances(symprec, angle_tolerance)

    matrix = validate_scaling_matrix(
        scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    )
    det = scaling_determinant(matrix)
    overrides = parse_set_occupancy(occupancy_overrides)
    data = read_cif(input_cif)
    data.atoms = _apply_overrides(data.atoms, overrides)
    data.atoms = _fill_and_validate_vacancies(data.atoms, vacancy_symbol)
    labels = [atom.label for atom in data.atoms]
    if len(set(labels)) != len(labels):
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        raise ValueError(
            "occupancy expansion produced duplicate CIF labels: "
            + ", ".join(duplicates)
        )

    expected_spacegroup = (
        parent_spacegroup if parent_spacegroup is not None
        else data.spacegroup_number
    )
    audit = audit_symmetry_grouping(
        data,
        expected_spacegroup=expected_spacegroup,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        symmetrize=symmetrize,
    )
    _check_integerizable_groups(data.atoms, audit["site_orbits"], det)
    _check_target_formula(
        data.atoms, audit["site_orbits"], det, vacancy_symbol, target_formula,
    )

    write_cif(data, output_cif)
    manifest_path = f"{output_cif}.manifest.json"
    orbit_path = f"{output_cif}.orbit_grouping.json"
    manifest = base_manifest(
        mode="strict_exhaustive" if strict else "shry_prepare",
        input_cif=os.path.abspath(output_cif),
        original_cif=os.path.abspath(input_cif),
        input_sha256=sha256_file(input_cif),
        prepared_sha256=sha256_file(output_cif),
        parent_spacegroup_declared=expected_spacegroup,
        parent_spacegroup_requested=parent_spacegroup,
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
    vacancy = element_symbol(vacancy_symbol)
    for group in groups.values():
        parsed = [(atom, parse_fraction(atom.occupancy)) for atom in group]
        for atom, occupancy in parsed:
            if element_symbol(atom.type_symbol) != vacancy:
                real_element_symbol(atom.type_symbol)
            if occupancy < 0 or occupancy > 1:
                raise ValueError(
                    f"occupancy for {atom.label} must be between 0 and 1, "
                    f"got {fraction_string(occupancy)}"
                )
        total = sum((occupancy for _, occupancy in parsed), Fraction(0))
        if total > 1:
            labels = ", ".join(atom.label for atom in group)
            raise ValueError(f"occupancy sum exceeds 1 for site {labels}: {total}")
        # A zero-occupancy row is not an atom. Keeping it can make downstream
        # disorder backends invent a species with zero sites, so omit it.
        out.extend(atom for atom, occupancy in parsed if occupancy > 0)
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
            raise ValueError(
                f"symmetry orbit {record['orbit_id']} representative {rep} "
                "is absent from the prepared atom rows"
            )
        occupancies: dict[str, Fraction] = {}
        for atom in group_atoms:
            species = element_symbol(atom.type_symbol)
            occupancies[species] = (
                occupancies.get(species, Fraction(0))
                + parse_fraction(atom.occupancy)
            )
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


def _check_target_formula(
    atoms: list[AtomRow],
    site_orbits: list[dict],
    determinant: int,
    vacancy_symbol: str,
    target_formula: str | None,
) -> None:
    """Verify the requested ordered-supercell formula before enumeration."""
    if target_formula is None:
        return
    expected = parse_formula(target_formula)
    atoms_by_coord: dict[tuple[str, str, str], list[AtomRow]] = {}
    for atom in atoms:
        atoms_by_coord.setdefault(atom.frac_key, []).append(atom)
    vacancy = element_symbol(vacancy_symbol)
    counts: dict[str, Fraction] = {}
    for record in site_orbits:
        representative = tuple(record["representative"])
        multiplier = int(record["multiplicity"]) * determinant
        for atom in atoms_by_coord.get(representative, []):
            species = element_symbol(atom.type_symbol)
            if species == vacancy:
                continue
            counts[species] = (
                counts.get(species, Fraction(0))
                + parse_fraction(atom.occupancy) * multiplier
            )
    non_integer = {species: value for species, value in counts.items()
                   if value.denominator != 1}
    if non_integer:
        raise ValueError(
            f"target-formula audit found non-integer ordered counts: {non_integer}"
        )
    actual = {species: value.numerator for species, value in counts.items()
              if value != 0}
    if actual != expected:
        raise ValueError(
            f"target formula mismatch before enumeration: got {actual}, "
            f"expected {expected} ({target_formula})"
        )
