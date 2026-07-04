"""Subprocess wrapper for the SHRY CLI."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass


@dataclass
class ShryResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str


class ShryBackend:
    """Run SHRY through its CLI, optionally from an isolated command path.

    SHRY pins ``pymatgen<=2023.10.4`` which conflicts with xtalkit's own
    enumeration extra, so SHRY lives in its own isolated environment (typically
    installed via ``uv tool install shry``). xtalkit reaches it only through
    this subprocess wrapper. The isolated interpreter is recovered from the
    ``shry`` executable's shebang so the SHRY-side version (§5) and pip freeze
    (§5/§8) reflect *that* environment, not xtalkit's.
    """

    def __init__(self, command: str | None = None):
        self.command = command or os.environ.get("XTALKIT_SHRY_CMD") or "shry"

    def resolve(self) -> str:
        exe = shutil.which(self.command)
        if exe is None:
            raise RuntimeError(
                "SHRY executable not found. Install SHRY in an isolated environment "
                "and set XTALKIT_SHRY_CMD to its CLI path, or make `shry` available "
                "on PATH."
            )
        return exe

    def python_path(self) -> str | None:
        """Return the SHRY env's python, parsed from the ``shry`` shebang."""
        try:
            exe = self.resolve()
        except RuntimeError:
            return None
        try:
            with open(exe, encoding="utf-8") as f:
                first = f.readline().strip()
        except (OSError, UnicodeDecodeError):
            return None
        if not first.startswith("#!"):
            return None
        interp = first[2:].strip().split()[0]
        return interp if interp and os.path.exists(interp) else None

    def run(self, args: list[str], cwd: str | None = None) -> ShryResult:
        exe = self.resolve()
        command = [exe, *args]
        result = subprocess.run(
            command,
            cwd=cwd,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "SHRY command failed:\n"
                f"  {' '.join(command)}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )
        return ShryResult(command, result.returncode, result.stdout, result.stderr)

    def version(self) -> str | None:
        """Return the SHRY package version from its isolated environment.

        SHRY's CLI has no working ``--version`` flag (argparse prints usage),
        so we ask the SHRY interpreter for ``shry.__version__`` / pip metadata.
        """
        py = self.python_path()
        if py:
            for probe in (
                ["-c", "import shry,sys;sys.stdout.write(shry.__version__)"],
                ["-c", "import importlib.metadata as m,sys;sys.stdout.write(m.version('shry'))"],
            ):
                try:
                    r = subprocess.run([py, *probe], capture_output=True,
                                       text=True, check=False)
                except OSError:
                    continue
                if r.returncode == 0 and r.stdout.strip():
                    return r.stdout.strip()
        # Fallback: try the CLI --version (older/custom builds).
        try:
            result = self.run(["--version"])
        except Exception:
            return None
        text = (result.stdout + "\n" + result.stderr).strip()
        return text or None

    def pip_freeze(self, path: str) -> bool:
        """Write the SHRY env's package snapshot to ``path``.

        Returns True if the SHRY interpreter was used, False if a fallback
        (xtalkit's own env) had to be used instead. Uses ``importlib.metadata``
        because uv tool envs (where SHRY lives) ship without ``pip``.
        """
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        from xtalkit.enumeration.manifest import _freeze_via
        py = self.python_path()
        if py:
            text = _freeze_via(py)
            if text and not text.startswith("pip freeze unavailable"):
                with open(path, "w", encoding="utf-8") as f:
                    f.write(text)
                return True
        # Fallback: the calling (xtalkit) interpreter.
        import sys
        with open(path, "w", encoding="utf-8") as f:
            f.write(_freeze_via(sys.executable))
        return False


def matrix_args(matrix: list[list[int]]) -> list[str]:
    """Flatten a 3x3 diagonal/general supercell matrix for SHRY CLI."""
    flat = [str(int(v)) for row in matrix for v in row]
    if matrix[0][1:] == [0, 0] and matrix[1][0] == 0 and matrix[1][2] == 0 \
            and matrix[2][:2] == [0, 0]:
        return [str(matrix[0][0]), str(matrix[1][1]), str(matrix[2][2])]
    return flat


def parse_count_output(stdout: str, stderr: str = "") -> int:
    """Extract SHRY count-only result from common CLI output forms.

    SHRY 1.1.x prints two numbers: the raw combination count
    ("Total number of combinations is N") and the symmetry-inequivalent count
    ("Expected unique patterns is M"). We want the latter — it is the number of
    structures `shry enum` will actually write. The earlier regexes are kept as
    fallbacks for other SHRY builds / mocked output.
    """
    text = stdout + "\n" + stderr
    patterns = [
        r"Expected\s+unique\s+patterns\s+is\s+(\d+)",
        r"unique\s+patterns\s+is\s+(\d+)",
        r"inequivalent(?:\s+structures?)?\s*[:=]\s*(\d+)",
        r"count(?:-only)?(?:\s+result)?\s*[:=]\s*(\d+)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return int(m.group(1))
    # Ambiguous fallback: only trust a lone integer (e.g. mocked "7").
    numbers = re.findall(r"\b\d+\b", text)
    if len(numbers) == 1:
        return int(numbers[0])
    raise ValueError(f"Could not parse SHRY count-only output:\n{text}")


def matrix_args(matrix: list[list[int]]) -> list[str]:
    """Flatten a 3x3 diagonal/general supercell matrix for SHRY CLI."""
    flat = [str(int(v)) for row in matrix for v in row]
    if matrix[0][1:] == [0, 0] and matrix[1][0] == 0 and matrix[1][2] == 0 \
            and matrix[2][:2] == [0, 0]:
        return [str(matrix[0][0]), str(matrix[1][1]), str(matrix[2][2])]
    return flat


def parse_count_output(stdout: str, stderr: str = "") -> int:
    """Extract SHRY count-only result from common CLI output forms.

    SHRY 1.1.x prints two numbers: the raw combination count
    ("Total number of combinations is N") and the symmetry-inequivalent count
    ("Expected unique patterns is M"). We want the latter — it is the number of
    structures `shry enum` will actually write. The earlier regexes are kept as
    fallbacks for other SHRY builds / mocked output.
    """
    text = stdout + "\n" + stderr
    patterns = [
        r"Expected\s+unique\s+patterns\s+is\s+(\d+)",
        r"unique\s+patterns\s+is\s+(\d+)",
        r"inequivalent(?:\s+structures?)?\s*[:=]\s*(\d+)",
        r"count(?:-only)?(?:\s+result)?\s*[:=]\s*(\d+)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return int(m.group(1))
    # Ambiguous fallback: only trust a lone integer (e.g. mocked "7").
    numbers = re.findall(r"\b\d+\b", text)
    if len(numbers) == 1:
        return int(numbers[0])
    raise ValueError(f"Could not parse SHRY count-only output:\n{text}")
