"""Environment setup for running enumlib through pymatgen.

``pymatgen.command_line.enumlib_caller`` shells out to the Fortran ``enum.x``
and ``makestr.x`` binaries, locating them with ``shutil.which`` at *module
import time*. So before that module is imported we must put the binaries on
PATH, and on Windows fix a few ``shutil.which`` / native-DLL quirks.

Binaries are found in the first existing location, in priority order:

1. ``$XTALKIT_ENUMLIB_BIN`` (a directory holding ``enum.x`` + ``makestr.x``)
2. a user-local default — ``~/.local/share/xtalkit/bin`` on Linux/macOS,
   ``%LOCALAPPDATA%/xtalkit/bin`` on Windows — which is where
   ``scripts/build_enumlib.sh`` installs them
3. an in-repo ``enumlib_src/enumlib/src`` clone (development; gitignored)

On Windows three extra workarounds are applied: ``PATHEXT`` is extended with
``.X``/``.PY``; the conda env's ``Library/bin`` DLL directories are registered
with ``os.add_dll_directory`` (Python 3.8+ requires this for native module
loading); and ``shutil.which`` is wrapped to return absolute paths only.

Caveat: the ``shutil.which`` absolute-path patch is global and persists for
the whole process (not only the enumerate path). That is acceptable for the
one-shot CLI; if you embed xtalkit as a library and depend on ``which``
returning relative paths, run ``enumerate_structures`` in a subprocess.
"""

from __future__ import annotations

import os
import shutil
import sys

_ALREADY_SETUP = False
_ORIG_WHICH = shutil.which


def _absolute_which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Wrap ``shutil.which`` to return absolute paths only."""
    result = _ORIG_WHICH(cmd, mode=mode, path=path)
    if result and not os.path.isabs(result):
        abs_result = os.path.abspath(result)
        if os.path.exists(abs_result):
            return abs_result
    return result


def _enumlib_bin_dirs() -> list[str]:
    """Candidate directories that may contain ``enum.x`` / ``makestr.x``.

    Ordered by priority; the first that exists and holds ``enum.x`` wins.
    """
    dirs: list[str] = []
    env_dir = os.environ.get("XTALKIT_ENUMLIB_BIN")
    if env_dir:
        dirs.append(os.path.abspath(os.path.expanduser(env_dir)))
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or os.path.join(
            os.path.expanduser("~"), "AppData", "Local")
        dirs.append(os.path.join(base, "xtalkit", "bin"))
    else:
        xdg = os.environ.get("XDG_DATA_HOME") or os.path.join(
            os.path.expanduser("~"), ".local", "share")
        dirs.append(os.path.join(xdg, "xtalkit", "bin"))
    # In-repo development clone (covered by .gitignore entry `enumlib_src/`).
    here = os.path.dirname(os.path.abspath(__file__))
    dirs.append(os.path.normpath(
        os.path.join(here, "..", "enumlib_src", "enumlib", "src")))
    return dirs


def _prepend_to_path(directory: str) -> None:
    """Prepend ``directory`` to PATH (no-op if already present)."""
    path = os.environ.get("PATH", "")
    parts = path.split(os.pathsep) if path else []
    if directory not in parts:
        os.environ["PATH"] = directory + os.pathsep + path


def setup_for_enumlib() -> None:
    """Make enumlib binaries discoverable + apply Windows workarounds.

    The binary-directory lookup runs every call (idempotent, cheap) so that
    binaries installed after the first call are still picked up. The
    Windows-only mutations and the ``shutil.which`` patch run once.

    Safe to call repeatedly.
    """
    # All platforms: put the compiled binaries on PATH before enumlib_caller's
    # import-time ``which("enum.x")`` runs.
    for d in _enumlib_bin_dirs():
        if os.path.isdir(d) and os.path.exists(os.path.join(d, "enum.x")):
            _prepend_to_path(d)
            break

    global _ALREADY_SETUP
    if _ALREADY_SETUP:
        return
    _ALREADY_SETUP = True

    if sys.platform == "win32":
        ext = os.environ.get("PATHEXT", "")
        for required in (".X", ".PY"):
            if required.upper() not in ext.upper():
                os.environ["PATHEXT"] = ext + (";" if ext else "") + required
                ext = os.environ["PATHEXT"]

        # When invoked directly (not via `conda activate`), the env's
        # Library/bin and Library/mingw-w64/bin may not be on PATH. Add them
        # so enum.x / makestr.x are discoverable and scipy's DLLs load.
        env_root = sys.prefix
        for sub in (os.path.join("Library", "bin"),
                    os.path.join("Library", "mingw-w64", "bin"),
                    os.path.join("Library", "usr", "bin"),
                    "Scripts"):
            d = os.path.join(env_root, sub)
            if os.path.isdir(d):
                try:
                    os.add_dll_directory(d)
                except (OSError, AttributeError):
                    pass
                _prepend_to_path(d)

    if shutil.which is not _absolute_which:
        shutil.which = _absolute_which
