"""Subprocess wrapper for the SHRY CLI."""

from __future__ import annotations

import os
import math
import re
import shlex
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
    installed via ``uv tool install --python 3.11 'shry==1.1.8'``). xtalkit
    reaches it only through this subprocess wrapper. The isolated interpreter
    is recovered from the ``shry`` executable's shebang so the SHRY-side
    version and package snapshot reflect *that* environment, not xtalkit's.
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
        shebang = shlex.split(first[2:].strip())
        if not shebang:
            return None
        interp = shebang[0]
        # Support portable ``#!/usr/bin/env python`` entry points as well as
        # uv's usual absolute-interpreter shebang.
        if os.path.basename(interp) == "env" and len(shebang) >= 2:
            interp = shutil.which(shebang[1]) or shebang[1]
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
                f"  {shlex.join(command)}\n"
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


def scaling_determinant(matrix: list[list[int]]) -> int:
    """Return the exact determinant of a validated 3x3 integer matrix."""
    a = matrix
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def validate_scaling_matrix(matrix: list[list[int]]) -> list[list[int]]:
    """Normalize a 3x3 integer supercell matrix and require positive volume."""
    try:
        rows = list(matrix)
        matrix_rows = [list(row) for row in rows]
    except TypeError as exc:
        raise ValueError("scaling matrix must be exactly 3x3") from exc
    if len(matrix_rows) != 3 or any(len(row) != 3 for row in matrix_rows):
        raise ValueError("scaling matrix must be exactly 3x3")
    normalized: list[list[int]] = []
    for row in matrix_rows:
        normalized_row: list[int] = []
        for value in row:
            if isinstance(value, bool):
                raise ValueError(f"scaling matrix entry {value!r} is not an integer")
            try:
                numeric = float(value)
                integer = int(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(
                    f"scaling matrix entry {value!r} is not a finite integer"
                ) from exc
            if not math.isfinite(numeric) or numeric != integer:
                raise ValueError(f"scaling matrix entry {value!r} is not an integer")
            normalized_row.append(integer)
        normalized.append(normalized_row)
    determinant = scaling_determinant(normalized)
    if determinant <= 0:
        raise ValueError(
            f"scaling matrix determinant must be positive, got {determinant}"
        )
    return normalized


def validate_tolerances(
    symprec: float,
    angle_tolerance: float,
    atol: float | None = None,
) -> None:
    """Validate finite tolerances shared by every SHRY workflow stage."""
    if not math.isfinite(symprec) or symprec <= 0:
        raise ValueError("symprec must be a finite positive number")
    if (not math.isfinite(angle_tolerance)
            or not (angle_tolerance == -1 or 0 <= angle_tolerance < 180)):
        raise ValueError(
            "angle_tolerance must be -1 (automatic) or a finite value in [0, 180)"
        )
    if atol is not None and (not math.isfinite(atol) or atol <= 0):
        raise ValueError("atol must be a finite positive number")


def matrix_args(matrix: list[list[int]]) -> list[str]:
    """Flatten a validated diagonal/general supercell matrix for SHRY CLI."""
    matrix = validate_scaling_matrix(matrix)
    flat = [str(v) for row in matrix for v in row]
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


def parse_total_combinations(stdout: str, stderr: str = "") -> int | None:
    """Return SHRY's raw combinatorial count when the backend reports it."""
    text = stdout + "\n" + stderr
    match = re.search(
        r"Total\s+number\s+of\s+combinations\s+is\s+(\d+)",
        text,
        flags=re.IGNORECASE,
    )
    return int(match.group(1)) if match else None
