"""Manifest helpers for audited SHRY workflows."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

from xtalkit import __version__


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: str, data: dict) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def append_jsonl(path: str, data: dict) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(data, sort_keys=True) + "\n")


def _dependency_versions() -> dict:
    """Best-effort version probe for the audit/dedup dependencies (plan §5).

    pymatgen (StructureMatcher dedup, §10) and spglib (symmetry audit, §3.4)
    are lazy-imported so the core toolkit still works without the ``enumerate``
    extra; missing packages record as ``null`` rather than failing the manifest.
    """
    versions: dict[str, str | None] = {}
    for name in ("pymatgen", "spglib"):
        try:
            mod = __import__(name)
        except ImportError:
            versions[name] = None
            continue
        # spglib exposes __version__; pymatgen does not (use importlib.metadata).
        ver = getattr(mod, "__version__", None)
        if ver is None:
            try:
                from importlib.metadata import version as _pkg_version
                ver = _pkg_version(name)
            except Exception:
                ver = None
        versions[name] = ver
    return versions


def base_manifest(**kwargs) -> dict:
    data = {
        "xtalkit_module": "shry",
        "mode": "strict_exhaustive",
        "backend": "SHRY",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "versions": {
            "xtalkit": __version__,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            **_dependency_versions(),
        },
        "sampling_enabled": False,
        "energy_filter_enabled": False,
    }
    data.update(kwargs)
    return data


def _freeze_via(python: str) -> str:
    """List installed packages with ``importlib.metadata`` (no pip needed).

    uv-managed environments (xtalkit's venv, SHRY's tool env) ship without
    ``pip``, so ``python -m pip freeze`` fails there. ``importlib.metadata`` is
    in the stdlib and works everywhere.
    """
    script = (
        "import importlib.metadata as m,sys;"
        "out=[(d.metadata['Name'] or '',d.version) for d in m.distributions()];"
        "sys.stdout.write('\\n'.join(f'{n}=={v}' for n,v in sorted(out)))"
    )
    try:
        r = subprocess.run([python, "-c", script],
                           capture_output=True, text=True, check=False)
    except OSError as exc:
        return f"pip freeze unavailable: {exc}\n"
    return r.stdout if r.returncode == 0 else f"pip freeze unavailable: {r.stderr}\n"


def pip_freeze(path: str) -> None:
    """Write a best-effort package snapshot of the current (xtalkit) env."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(_freeze_via(sys.executable))
