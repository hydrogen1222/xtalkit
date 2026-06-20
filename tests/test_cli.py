import os
import subprocess
import sys
import tempfile

import pytest


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


def test_cli_skeleton_no_output_flag(tmp_path):
    """Without -o, output goes to current directory with default name."""
    import os
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        result = run_xtalkit(
            "skeleton", "--sg", "216", "--wyckoff", "4a",
            "--format", "cif",
        )
        # Should succeed and create file
        assert result.returncode == 0
        files = list(tmp_path.glob("*.cif"))
        assert len(files) > 0
    finally:
        os.chdir(cwd)


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


def test_cli_no_args_enters_tui():
    """With no args, warn that TUI is not implemented yet (Task 8)."""
    result = run_xtalkit()
    # Should not crash
    assert result.returncode == 0
    output = result.stdout + result.stderr
    assert "xtalkit" in output.lower() or "menu" in output.lower() or "tui" in output.lower()
