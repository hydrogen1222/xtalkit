import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest

from xtalkit.cli import _path_arg


@pytest.fixture
def simple_cif():
    return os.path.join(os.path.dirname(__file__), "fixtures", "simple.cif")


def run_xtalkit(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "xtalkit.cli", *args],
        capture_output=True,
        text=True,
    )


def test_cli_version():
    result = run_xtalkit("--version")
    assert result.returncode == 0
    assert "xtalkit" in result.stdout


def test_cli_top_level_help_points_to_a_real_first_step():
    result = run_xtalkit("--help")
    assert result.returncode == 0
    assert "xtalkit info --sg 216" in result.stdout
    assert '--atom "Na 4a"' in result.stdout
    assert "without a subcommand" in result.stdout
    assert "README.zh-CN.md" in result.stdout


def test_cli_info():
    result = run_xtalkit("info", "--sg", "216")
    assert result.returncode == 0
    assert "216" in result.stdout
    assert "Wyckoff" in result.stdout or "4a" in result.stdout


def test_cli_info_invalid_sg():
    result = run_xtalkit("info", "--sg", "999")
    assert result.returncode != 0


def test_cli_skeleton():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "skel")
        result = run_xtalkit(
            "skeleton", "--sg", "216", "--wyckoff", "4a,4c",
            "--format", "cif", "-o", out,
        )
        assert result.returncode == 0
        assert os.path.exists(out + ".cif")


def test_cli_skeleton_no_output_flag(tmp_path, monkeypatch):
    """Without -o, output goes to current directory with default name."""
    monkeypatch.chdir(tmp_path)
    result = run_xtalkit(
        "skeleton", "--sg", "216", "--wyckoff", "4a",
        "--format", "cif",
    )
    # Should succeed and create file
    assert result.returncode == 0
    files = list(tmp_path.glob("*.cif"))
    assert len(files) > 0


def test_cli_mark_help():
    result = run_xtalkit("mark", "--help")
    assert result.returncode == 0
    assert "--sg" in result.stdout
    assert "--wyckoff" in result.stdout
    assert "--mode" in result.stdout


def test_cli_skeleton_help():
    result = run_xtalkit("skeleton", "--help")
    assert result.returncode == 0
    assert "--sg" in result.stdout


def test_cli_info_help():
    result = run_xtalkit("info", "--help")
    assert result.returncode == 0


def test_cli_no_args_in_noninteractive_shell_prints_help():
    """Redirected/CI invocation must not block waiting for TUI input."""
    result = run_xtalkit()
    assert result.returncode == 0
    assert "usage: xtalkit" in result.stdout
    assert "Typical first steps" in result.stdout


def test_cli_path_argument_expands_home_directory():
    assert _path_arg("~/example.cif") == str(Path.home() / "example.cif")
    assert _path_arg("relative/example.cif") == "relative/example.cif"


def test_cli_fetch():
    """'fetch' subcommand verifies supported space group data is intact.

    All 230 SGs are populated from the bundled dataset; fetch must succeed
    (exit 0) and report 230/230.
    """
    result = run_xtalkit("fetch")
    assert result.returncode == 0
    assert "supported" in result.stdout.lower()
    assert "230/230" in result.stdout


def test_cli_build_nacl():
    """Build NaCl (Fm-3m) from refinement parameters via flags."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "NaCl")
        result = run_xtalkit(
            "build", "--sg", "225",
            "--cell", "5.64 5.64 5.64 90 90 90",
            "--atom", "Na 4a", "--atom", "Cl 4b",
            "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "NaCl" in result.stdout
        assert os.path.exists(out + ".cif")


def test_cli_build_free_param_and_xyz():
    """Build with a free coordinate + xyz output."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "li")
        result = run_xtalkit(
            "build", "--sg", "216",
            "--cell", "5.9 5.9 5.9 90 90 90",
            "--atom", "Li 16e 0.3",
            "--format", "cif,xyz", "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert os.path.exists(out + ".cif")
        assert os.path.exists(out + ".xyz")


def test_cli_build_partial_occupancy():
    """Partial occupancy is accepted and reported."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "dis")
        result = run_xtalkit(
            "build", "--sg", "225", "--cell", "4 4 4 90 90 90",
            "--atom", "Li 4a 0.5", "--atom", "Cu 4a 0.5",
            "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "LiCu" in result.stdout


def test_cli_build_spec_json(tmp_path):
    """Build from a JSON spec file."""
    spec = tmp_path / "spec.json"
    spec.write_text(
        '{"sg":225,"cell":{"a":5.64,"b":5.64,"c":5.64,'
        '"alpha":90,"beta":90,"gamma":90},'
        '"atoms":[{"element":"Na","wyckoff":"4a","occ":1.0},'
        '{"element":"Cl","wyckoff":"4b","occ":1.0}]}'
    )
    out = str(tmp_path / "fromspec")
    result = run_xtalkit("build", "--spec", str(spec), "-o", out)
    assert result.returncode == 0, result.stderr
    assert os.path.exists(out + ".cif")


def test_cli_build_bad_free_count():
    """Wrong number of free values is a clean error (exit 1)."""
    result = run_xtalkit(
        "build", "--sg", "216", "--cell", "5 5 5 90 90 90",
        "--atom", "Li 16e 0.3 0.4 0.5",
    )
    assert result.returncode != 0
    assert "expects" in result.stderr


def test_cli_build_help():
    result = run_xtalkit("build", "--help")
    assert result.returncode == 0
    assert "--sg" in result.stdout
    assert "--atom" in result.stdout
    assert "--spec" in result.stdout


def test_cli_ewald_help():
    result = run_xtalkit("ewald", "--help")
    assert result.returncode == 0
    assert "--layout" in result.stdout
    assert "--top-n" in result.stdout
    assert "--group" in result.stdout
    assert "--jobs" in result.stdout


def test_cli_build_atom_frac_mode():
    """Build from direct fractional coordinates (--atom-frac)."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "NaCl")
        result = run_xtalkit(
            "build", "--sg", "225", "--cell", "5.64 5.64 5.64 90 90 90",
            "--atom-frac", "Na 0 0 0", "--atom-frac", "Cl 0.5 0.5 0.5",
            "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "NaCl" in result.stdout
        # detect_wyckoff reports the orbit label
        assert "4a" in result.stdout and "4b" in result.stdout
        assert os.path.exists(out + ".cif")


def test_cli_build_atom_frac_noncanonical_rep():
    """Non-canonical fractional coords are accepted (SG 137 4d via 4-fold)."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "li")
        result = run_xtalkit(
            "build", "--sg", "137", "--cell", "8.69 8.69 12.6 90 90 90",
            "--atom-frac", "Li 0 0.5 0.9446 1", "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "4d" in result.stdout  # detected despite non-canonical (0,1/2,z)


def test_cli_build_atom_frac_skips_zero_occupancy():
    """occ=0 atoms (absent) are skipped with a note, not an error."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "x")
        result = run_xtalkit(
            "build", "--sg", "225", "--cell", "4 4 4 90 90 90",
            "--atom-frac", "Li 0 0 0 0.5", "--atom-frac", "Cu 0 0 0 0",
            "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "skipping" in result.stderr


# Li6PS5Cl-like argyrodite (F-43m) with Cl/S mixing on 4a and 4d and
# half-occupied Li on 48h: Li24 P4 S20 Cl4 per cell -> Li6PS5Cl.
_ARGYRODITE_SPEC = {
    "sg": 216,
    "cell": {"a": 9.86, "b": 9.86, "c": 9.86,
             "alpha": 90, "beta": 90, "gamma": 90},
    "atoms": [
        {"element": "Li", "xyz": [0.182, 0.182, 0.022], "occ": 0.5},
        {"element": "P", "xyz": [0.5, 0.5, 0.5]},
        {"element": "S", "wyckoff": "16e", "free": [0.62]},
        {"element": "Cl", "xyz": [0, 0, 0], "occ": 0.5},
        {"element": "S", "xyz": [0, 0, 0], "occ": 0.5},
        {"element": "S", "xyz": [0.75, 0.75, 0.75], "occ": 0.5},
        {"element": "Cl", "xyz": [0.75, 0.75, 0.75], "occ": 0.5},
    ],
}


def test_cli_build_spec_with_fractional_coords(tmp_path):
    """A JSON spec may give fractional 'xyz' and mix it with Wyckoff rows."""
    import json
    spec = tmp_path / "lpsc.json"
    spec.write_text(json.dumps(_ARGYRODITE_SPEC))
    out = str(tmp_path / "lpsc")
    result = run_xtalkit("build", "--spec", str(spec), "-o", out)
    assert result.returncode == 0, result.stderr
    assert "Li6PS5Cl" in result.stdout
    assert "48h" in result.stdout and "16e" in result.stdout
    assert os.path.exists(out + ".cif")


def test_cli_build_spec_rejects_xyz_and_wyckoff_together(tmp_path):
    import json
    spec = tmp_path / "bad.json"
    spec.write_text(json.dumps({
        "sg": 225,
        "cell": {"a": 5.64, "b": 5.64, "c": 5.64,
                 "alpha": 90, "beta": 90, "gamma": 90},
        "atoms": [{"element": "Na", "wyckoff": "4a", "xyz": [0, 0, 0]}],
    }))
    result = run_xtalkit("build", "--spec", str(spec))
    assert result.returncode != 0
    assert "not both" in result.stderr


def test_cli_build_mixes_atom_and_atom_frac():
    """--atom rows are kept (not silently dropped) next to --atom-frac."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "NaCl")
        result = run_xtalkit(
            "build", "--sg", "225", "--cell", "5.64 5.64 5.64 90 90 90",
            "--atom", "Na 4a", "--atom-frac", "Cl 0.5 0.5 0.5",
            "-o", out,
        )
        assert result.returncode == 0, result.stderr
        assert "NaCl" in result.stdout
        assert "4a" in result.stdout and "4b" in result.stdout


def test_cli_build_overoccupied_site_hints_at_rounding():
    result = run_xtalkit(
        "build", "--sg", "225", "--cell", "4 4 4 90 90 90",
        "--atom-frac", "Cl 0 0 0 0.616", "--atom-frac", "S 0 0 0 0.385",
    )
    assert result.returncode != 0
    assert "renormalize" in result.stderr
