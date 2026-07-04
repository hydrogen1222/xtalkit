"""Subprocess wrapper for the independent ``supercell`` backend."""

from __future__ import annotations

import os
import re
import shutil
import subprocess


def run_supercell_count(
    cif_path: str,
    scaling_matrix: list[list[int]] | None = None,
    command: str | None = None,
) -> int:
    """Run the external supercell program and parse an inequivalent count.

    The exact supercell CLI differs across builds, so the command path can be
    configured with ``XTALKIT_SUPERCELL_CMD``. This wrapper is intentionally
    conservative: it only reports a count if the output has an unambiguous
    integer.
    """
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)
    exe_name = command or os.environ.get("XTALKIT_SUPERCELL_CMD") or "supercell"
    exe = shutil.which(exe_name)
    if exe is None:
        raise RuntimeError(
            "supercell executable not found. Install the independent supercell "
            "program and set XTALKIT_SUPERCELL_CMD, or put `supercell` on PATH."
        )

    matrix = scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    # Common supercell builds accept input file plus -s/--supercell-like sizes.
    # Users with different builds can wrap the executable path in a small script.
    args = [exe, cif_path, "--supercell", *[str(v) for row in matrix for v in row], "--count-only"]
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "supercell command failed:\n"
            f"  {' '.join(args)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return parse_supercell_count(result.stdout + "\n" + result.stderr)


def parse_supercell_count(output: str) -> int:
    """Parse a count from supercell output."""
    patterns = [
        r"inequivalent(?:\s+structures?)?\s*[:=]\s*(\d+)",
        r"unique(?:\s+structures?)?\s*[:=]\s*(\d+)",
        r"total(?:\s+structures?)?\s*[:=]\s*(\d+)",
        r"count(?:\s+result)?\s*[:=]\s*(\d+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, output, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    numbers = re.findall(r"\b\d+\b", output)
    if len(numbers) == 1:
        return int(numbers[0])
    raise ValueError(f"Could not parse supercell count output:\n{output}")
