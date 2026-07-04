import csv
from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace

from xtalkit.ewald import (
    EwaldRow,
    batch_ewald,
    group_rows,
    iter_structure_entries,
    parse_charges,
    split_rows,
    write_csv,
)


def test_parse_charges():
    charges = parse_charges(["Li:1", "Ge:4", "P:-5", "S:-2"])
    assert charges == {"Li": 1.0, "Ge": 4.0, "P": -5.0, "S": -2.0}


def test_iter_structure_entries_flat_and_nested(tmp_path):
    flat = tmp_path / "flat"
    flat.mkdir()
    (flat / "1.cif").write_text("data_1")
    (flat / "2.cif").write_text("data_2")
    (flat / "ignore.txt").write_text("nope")
    (flat / "nested").mkdir()
    (flat / "nested" / "3.cif").write_text("data_3")

    flat_entries = list(iter_structure_entries([str(flat)], layout="flat"))
    assert [rel for _, _, rel in flat_entries] == ["1.cif", "2.cif"]

    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "root.cif").write_text("root")
    sub = nested / "1"
    sub.mkdir()
    (sub / "1.cif").write_text("sub1")
    (sub / "2.cif").write_text("sub2")
    (sub / "skip.txt").write_text("skip")

    nested_entries = list(iter_structure_entries([str(nested)], layout="nested"))
    assert [rel for _, _, rel in nested_entries] == ["root.cif", "1/1.cif", "1/2.cif"]


def test_split_and_group_rows(tmp_path):
    src_a = tmp_path / "src_a"
    src_b = tmp_path / "src_b"
    src_a.mkdir()
    src_b.mkdir()
    (src_a / "1.cif").write_text("a1")
    (src_a / "2.cif").write_text("a2")
    (src_b / "1.cif").write_text("b1")

    rows = [
        EwaldRow(str(src_a), "1.cif", str(src_a / "1.cif"), "A", 1, 1.0),
        EwaldRow(str(src_b), "1.cif", str(src_b / "1.cif"), "B", 1, 2.0),
        EwaldRow(str(src_a), "2.cif", str(src_a / "2.cif"), "C", 1, 3.0),
    ]

    selected, remaining = split_rows(rows, top_n=2)
    assert [row.ewald_energy for row in selected] == [1.0, 2.0]
    assert [row.ewald_energy for row in remaining] == [3.0]

    info = group_rows(selected, remaining, str(tmp_path / "ewald_out"))
    selected_dir = Path(info["selected_dir"])
    remaining_dir = Path(info["remaining_dir"])

    assert (selected_dir / "src_a" / "1.cif").exists()
    assert (selected_dir / "src_b" / "1.cif").exists()
    assert (remaining_dir / "src_a" / "2.cif").exists()


def test_write_csv(tmp_path):
    rows = [
        EwaldRow("root", "1.cif", "/tmp/1.cif", "A", 2, -1.23),
        EwaldRow("root", "2.cif", "/tmp/2.cif", "B", 3, -0.5),
    ]
    out = tmp_path / "ranking.csv"
    write_csv(rows, str(out))

    with out.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        data = list(reader)

    assert data[0]["rank"] == "1"
    assert data[0]["relative_path"] == "1.cif"
    assert data[1]["ewald_energy"] == "-0.5"


def test_batch_ewald_parallel_branch(monkeypatch):
    entries = [("root", "/tmp/a.cif", "a.cif"), ("root", "/tmp/b.cif", "b.cif")]
    seen = {}

    class FakeExecutor:
        def __init__(self, max_workers):
            seen["max_workers"] = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, *args, **kwargs):
            fut = Future()
            fut.set_result(fn(*args, **kwargs))
            return fut

    def fake_score(source_root, path, relative_path, charges, guess, per_atom):
        energy = 2.0 if path.endswith("a.cif") else 1.0
        return EwaldRow(source_root, relative_path, path, "X", 1, energy)

    monkeypatch.setattr("xtalkit.ewald.iter_structure_entries", lambda paths, layout="flat": entries)
    monkeypatch.setattr("xtalkit.ewald._score_structure_entry", fake_score)
    monkeypatch.setattr("xtalkit.ewald.cf.ProcessPoolExecutor", FakeExecutor)

    rows = batch_ewald(["/tmp"], jobs=4)

    assert seen["max_workers"] == 2
    assert [row.path for row in rows] == ["/tmp/b.cif", "/tmp/a.cif"]
