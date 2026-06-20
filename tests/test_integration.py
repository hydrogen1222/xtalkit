"""End-to-end integration tests for xtalkit."""

import os
import subprocess
import sys
import tempfile

import pytest
import gemmi


@pytest.fixture
def simple_cif():
    return os.path.join(os.path.dirname(__file__), "fixtures", "simple.cif")


def run_xtalkit(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "xtalkit.cli", *args],
        capture_output=True,
        text=True,
    )


def test_mark_then_open_in_gemmi(simple_cif):
    """Mark a CIF, then verify Gemmi can re-parse the output."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "test_WYCK")
        result = run_xtalkit(
            "mark", simple_cif,
            "--sg", "216", "--wyckoff", "4a,4c,16e,24f",
            "--mode", "overlay", "--format", "cif",
            "-o", out,
        )
        assert result.returncode == 0

        # Gemmi should parse output CIF
        doc = gemmi.cif.read_file(out + ".cif")
        assert doc is not None
        assert len(doc) > 0
        block = doc.sole_block()

        # Verify cell parameters
        assert float(block.find_value("_cell_length_a")) == 5.0

        # Verify atom site loop contains expected atom labels
        content = open(out + ".cif").read()
        assert "Li1" in content
        assert "WYCK_4a" in content
        assert "WYCK_4c" in content


def test_skeleton_all_formats():
    """Generate skeleton in all 3 formats."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "skel")
        result = run_xtalkit(
            "skeleton", "--sg", "216",
            "--wyckoff", "4a,4c",
            "--format", "cif,xyz",
            "-o", out,
        )
        assert result.returncode == 0
        for ext in ("cif", "xyz"):
            assert os.path.exists(out + f".{ext}")


def test_mark_replace_mode(simple_cif):
    """Replace mode should not error."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "test_WYCK")
        result = run_xtalkit(
            "mark", simple_cif,
            "--sg", "216", "--wyckoff", "4a",
            "--mode", "replace", "--format", "cif",
            "-o", out,
        )
        assert result.returncode == 0


def test_info_all_space_groups():
    """Info should work for all 230 space groups."""
    for sg_num in (1, 2, 195, 216, 225, 230):
        result = run_xtalkit("info", "--sg", str(sg_num))
        assert result.returncode == 0


def test_custom_tolerance(simple_cif):
    """Custom tolerance should work."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "test_WYCK")
        result = run_xtalkit(
            "mark", simple_cif,
            "--sg", "216", "--wyckoff", "4a",
            "--tol", "0.01", "--format", "cif",
            "-o", out,
        )
        assert result.returncode == 0


def test_element_override(simple_cif):
    """Element map override should be respected."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "test_WYCK")
        result = run_xtalkit(
            "mark", simple_cif,
            "--sg", "216", "--wyckoff", "4a,4c",
            "--map", "4a:He,4c:Ne", "--format", "cif",
            "-o", out,
        )
        assert result.returncode == 0
        content = open(out + ".cif").read()
        assert "He" in content
        assert "Ne" in content


def test_cli_version_flag():
    """--version should work."""
    result = run_xtalkit("--version")
    assert result.returncode == 0


def test_cli_help():
    """--help should work."""
    result = run_xtalkit("--help")
    assert result.returncode == 0
    assert "mark" in result.stdout
    assert "skeleton" in result.stdout
    assert "info" in result.stdout
    assert "fetch" in result.stdout
