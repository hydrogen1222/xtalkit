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
        },
        "sampling_enabled": False,
        "energy_filter_enabled": False,
    }
    data.update(kwargs)
    return data


def pip_freeze(path: str) -> None:
    """Write a best-effort pip freeze snapshot."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            check=False,
        )
        text = result.stdout if result.returncode == 0 else result.stderr
    except Exception as exc:
        text = f"pip freeze unavailable: {exc}\n"
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
