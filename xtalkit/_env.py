"""Windows-specific environment fixes for pymatgen + enumlib.

Three issues need workarounds when running enumlib through pymatgen on Windows:

1. ``shutil.which`` ignores files whose extension is not in ``PATHEXT``.
   ``enum.x`` and ``makestr.x`` have a ``.x`` extension that isn't in the
   default PATHEXT, so we append ``.X`` (and ``.PY`` for makeStr.py fallback).

2. ``scipy.signal`` (transitively imported by pymatgen via
   ``pymatgen.electronic_structure.dos``) loads a native extension that
   depends on DLLs in ``Library/bin`` of the conda env. When Python is
   launched via ``conda run`` the env's ``Library/bin`` is on PATH but
   not registered with ``os.add_dll_directory`` (which Python 3.8+
   requires on Windows for native module loading).

3. ``shutil.which`` returns the first match, which on Windows can be a
   relative path (``.\\makestr.x``) when the current directory contains
   the file. ``subprocess.Popen`` cannot launch relative paths without
   ``shell=True``. We patch ``which`` to return only absolute paths.

The fixes are idempotent and only affect the ``enumerate`` code path.
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


def setup_for_enumlib() -> None:
    """Apply the three Windows workarounds. Safe to call repeatedly."""
    global _ALREADY_SETUP
    if _ALREADY_SETUP:
        return

    if sys.platform == "win32":
        ext = os.environ.get("PATHEXT", "")
        for required in (".X", ".PY"):
            if required.upper() not in ext.upper():
                os.environ["PATHEXT"] = ext + (";" if ext else "") + required
                ext = os.environ["PATHEXT"]

        # When invoked directly (not via `conda activate`), the env's
        # Library/bin and Library/mingw-w64/bin may not be on PATH.
        # Add them so enum.x / makestr.x are discoverable.
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
                path_list = os.environ.get("PATH", "").split(os.pathsep)
                if d not in path_list:
                    os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")

    if shutil.which is not _absolute_which:
        shutil.which = _absolute_which

    _ALREADY_SETUP = True
