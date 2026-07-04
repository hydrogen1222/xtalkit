"""Batch Ewald electrostatic energy ranking and grouping."""

from __future__ import annotations

import csv
import concurrent.futures as cf
import os
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class EwaldRow:
    """One ranked structure entry."""

    source_root: str
    relative_path: str
    path: str
    formula: str
    n_atoms: int
    ewald_energy: float

    def as_dict(self) -> dict[str, object]:
        return {
            "source_root": self.source_root,
            "relative_path": self.relative_path,
            "path": self.path,
            "formula": self.formula,
            "n_atoms": self.n_atoms,
            "ewald_energy": self.ewald_energy,
        }


def parse_charges(tokens: list[str]) -> dict[str, float]:
    """Parse ``Li:1 Ge:4 P:-5 S:-2`` into a charge map."""
    charges: dict[str, float] = {}
    for token in tokens:
        if ":" not in token:
            raise ValueError(
                f"bad charge token {token!r}; expected Element:Oxidation, e.g. Li:1"
            )
        el, val = token.split(":", 1)
        el = el.strip()
        if not el:
            raise ValueError(f"bad charge token {token!r}: empty element")
        try:
            charges[el] = float(val)
        except ValueError as exc:
            raise ValueError(f"bad charge token {token!r}: {val!r} is not a number") from exc
    return charges


def compute_ewald_energy(
    structure,
    charges: dict[str, float] | None = None,
    guess: bool = False,
    per_atom: bool = False,
) -> float:
    """Compute the Ewald electrostatic energy for one pymatgen structure."""
    try:
        from pymatgen.analysis.ewald import EwaldSummation
        from pymatgen.core.periodic_table import Species
    except ImportError as exc:  # pragma: no cover - exercised in integration envs
        raise RuntimeError(
            "pymatgen is required for Ewald scoring; run 'uv sync --extra enumerate'"
        ) from exc

    s = structure.copy()
    _assign_oxidation_states(s, charges=charges, guess=guess, specie_type=Species)
    energy = float(EwaldSummation(s).total_energy)
    return energy / len(s) if per_atom else energy


def _assign_oxidation_states(structure, charges, guess: bool, specie_type) -> None:
    """Attach oxidation states to a structure in place."""
    from pymatgen.core import DummySpecie

    if any(isinstance(sp, DummySpecie) for sp in structure.species):
        raise ValueError(
            "structure contains a vacancy/dummy species (X); Ewald needs ordered "
            "structures without vacancy pseudo-species"
        )

    if charges:
        present = {
            getattr(el, "symbol", str(el))
            for el in structure.composition.elements
        }
        missing = sorted(present - set(charges))
        if missing:
            raise ValueError(f"missing --charges for species: {missing}")
        structure.add_oxidation_state_by_element(charges)
        return

    if guess:
        structure.add_oxidation_state_by_guess()
        return

    if not any(isinstance(sp, specie_type) for sp in structure.species):
        raise ValueError(
            "structure has no oxidation states; pass --charges EL:OX ... or use --guess"
        )


def _is_structure_file(name: str) -> bool:
    """Return True for file names that pymatgen commonly reads as structures."""
    upper = name.upper()
    return name.lower().endswith(".cif") or upper.startswith("POSCAR") or upper.startswith("CONTCAR")


def _iter_flat_entries(paths: list[str]):
    """Yield files directly under each directory, or any explicit file path."""
    seen: set[str] = set()
    for raw in paths:
        path = Path(raw).expanduser()
        if not path.exists():
            raise FileNotFoundError(raw)
        if path.is_file():
            resolved = str(path.resolve())
            if resolved not in seen:
                seen.add(resolved)
                yield str(path.parent.resolve()), str(path), path.name
            continue

        for child in sorted(path.iterdir()):
            if child.is_file() and _is_structure_file(child.name):
                resolved = str(child.resolve())
                if resolved in seen:
                    continue
                seen.add(resolved)
                yield str(path.resolve()), str(child), child.name


def _iter_nested_entries(paths: list[str]):
    """Yield files one directory below each root directory.

    The root directory itself is also scanned for direct structure files.
    """
    seen: set[str] = set()
    for raw in paths:
        path = Path(raw).expanduser()
        if not path.exists():
            raise FileNotFoundError(raw)
        if path.is_file():
            resolved = str(path.resolve())
            if resolved not in seen:
                seen.add(resolved)
                yield str(path.parent.resolve()), str(path), path.name
            continue

        root = str(path.resolve())
        children = sorted(path.iterdir())
        for child in children:
            if child.is_file() and _is_structure_file(child.name):
                resolved = str(child.resolve())
                if resolved in seen:
                    continue
                seen.add(resolved)
                yield root, str(child), child.name
        for child in children:
            if not child.is_dir():
                continue
            for grandchild in sorted(child.iterdir()):
                if not grandchild.is_file() or not _is_structure_file(grandchild.name):
                    continue
                resolved = str(grandchild.resolve())
                if resolved in seen:
                    continue
                seen.add(resolved)
                yield root, str(grandchild), str(Path(child.name) / grandchild.name)


def iter_structure_entries(paths: list[str], layout: str = "flat"):
    """Yield ``(source_root, path, relative_path)`` triples for the batch run."""
    layout = layout.lower()
    if layout == "flat":
        yield from _iter_flat_entries(paths)
        return
    if layout == "nested":
        yield from _iter_nested_entries(paths)
        return
    raise ValueError("--layout must be 'flat' or 'nested'")


def batch_ewald(
    paths: list[str],
    charges: dict[str, float] | None = None,
    guess: bool = False,
    per_atom: bool = False,
    layout: str = "flat",
    jobs: int = 1,
) -> list[EwaldRow]:
    """Compute Ewald energy for every structure under ``paths`` and rank them."""
    try:
        from pymatgen.core import Structure
    except ImportError as exc:  # pragma: no cover - exercised in integration envs
        raise RuntimeError(
            "pymatgen is required for Ewald scoring; run 'uv sync --extra enumerate'"
        ) from exc

    rows: list[EwaldRow] = []
    entries = list(iter_structure_entries(paths, layout=layout))
    if not entries:
        raise ValueError(
            f"no structure files found under {paths} with layout={layout!r} "
            "(looked for CIF/POSCAR/CONTCAR files)"
        )

    if jobs < 0:
        raise ValueError("--jobs must be >= 0")
    workers = (os.cpu_count() or 1) if jobs == 0 else jobs
    workers = max(1, min(workers, len(entries)))
    if workers == 1 or len(entries) == 1:
        for source_root, path, relative_path in entries:
            rows.append(_score_structure_entry(
                source_root,
                path,
                relative_path,
                charges=charges,
                guess=guess,
                per_atom=per_atom,
            ))
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    _score_structure_entry,
                    source_root,
                    path,
                    relative_path,
                    charges,
                    guess,
                    per_atom,
                )
                for source_root, path, relative_path in entries
            ]
            for future in cf.as_completed(futures):
                try:
                    rows.append(future.result())
                except Exception as exc:  # noqa: BLE001
                    raise RuntimeError(str(exc)) from exc

    rows.sort(key=lambda row: row.ewald_energy)
    return rows


def _score_structure_entry(
    source_root: str,
    path: str,
    relative_path: str,
    charges: dict[str, float] | None,
    guess: bool,
    per_atom: bool,
) -> EwaldRow:
    """Score one structure file, suitable for multiprocessing."""
    try:
        from pymatgen.core import Structure
    except ImportError as exc:  # pragma: no cover - exercised in integration envs
        raise RuntimeError(
            "pymatgen is required for Ewald scoring; run 'uv sync --extra enumerate'"
        ) from exc

    try:
        struct = Structure.from_file(path)
        energy = compute_ewald_energy(
            struct,
            charges=charges,
            guess=guess,
            per_atom=per_atom,
        )
        return EwaldRow(
            source_root=source_root,
            relative_path=relative_path,
            path=path,
            formula=str(struct.composition.reduced_formula),
            n_atoms=len(struct),
            ewald_energy=energy,
        )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"{path}: {exc}") from exc


def sort_rows(rows: list[EwaldRow], descending: bool = False) -> list[EwaldRow]:
    """Return rows sorted by Ewald energy."""
    return sorted(rows, key=lambda row: row.ewald_energy, reverse=descending)


def split_rows(rows: list[EwaldRow], top_n: int | None = None) -> tuple[list[EwaldRow], list[EwaldRow]]:
    """Split rows into selected and remaining groups after sorting."""
    if top_n is None:
        return list(rows), []
    if top_n <= 0:
        raise ValueError("--top-n must be a positive integer")
    return list(rows[:top_n]), list(rows[top_n:])


def _source_root_labels(rows: list[EwaldRow]) -> dict[str, str]:
    """Build stable folder labels for one or more source roots."""
    roots = list(dict.fromkeys(row.source_root for row in rows))
    labels: dict[str, str] = {}
    seen: dict[str, int] = {}
    for root in roots:
        base = Path(root).name or "root"
        seen[base] = seen.get(base, 0) + 1
        labels[root] = base if seen[base] == 1 else f"{base}_{seen[base]}"
    return labels


def _copy_or_move(src: Path, dst: Path, move: bool) -> None:
    """Copy or move a file to ``dst``, creating parents first."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dst.resolve():
        return
    if move:
        shutil.move(str(src), str(dst))
    else:
        shutil.copy2(str(src), str(dst))


def group_rows(
    selected: list[EwaldRow],
    remaining: list[EwaldRow],
    output_dir: str,
    selected_name: str = "selected",
    remaining_name: str = "rest",
    move: bool = False,
) -> dict[str, str]:
    """Copy or move ranked structures into two folders."""
    out = Path(output_dir)
    selected_dir = out / selected_name
    remaining_dir = out / remaining_name
    selected_dir.mkdir(parents=True, exist_ok=True)
    remaining_dir.mkdir(parents=True, exist_ok=True)

    combined = selected + remaining
    if not combined:
        return {
            "selected_dir": str(selected_dir),
            "remaining_dir": str(remaining_dir),
        }

    labels = _source_root_labels(combined)
    use_root_labels = len(labels) > 1

    def _dest_base(row: EwaldRow, root_dir: Path) -> Path:
        rel = Path(row.relative_path)
        if use_root_labels:
            rel = Path(labels[row.source_root]) / rel
        return root_dir / rel

    for row in selected:
        _copy_or_move(Path(row.path), _dest_base(row, selected_dir), move)
    for row in remaining:
        _copy_or_move(Path(row.path), _dest_base(row, remaining_dir), move)

    return {
        "selected_dir": str(selected_dir),
        "remaining_dir": str(remaining_dir),
    }


def write_csv(rows: list[EwaldRow], path: str) -> None:
    """Write ranked rows to a CSV file."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "rank",
                "source_root",
                "relative_path",
                "path",
                "formula",
                "n_atoms",
                "ewald_energy",
            ],
        )
        writer.writeheader()
        for i, row in enumerate(rows, 1):
            d = row.as_dict()
            d["rank"] = i
            writer.writerow(d)
