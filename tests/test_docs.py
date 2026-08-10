"""Documentation contracts for the human-facing manuals."""

from pathlib import Path
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]
MANUALS = [
    ROOT / "README.md",
    ROOT / "README.zh-CN.md",
    ROOT / "docs/user/shry-enumeration.md",
    ROOT / "docs/user/shry-enumeration.zh-CN.md",
]


@pytest.mark.parametrize("manual", MANUALS)
def test_manual_local_links_exist(manual):
    text = manual.read_text(encoding="utf-8")
    for raw_target in re.findall(r"\[[^]]*\]\(([^)]+)\)", text):
        target = raw_target.strip().strip("<>").split("#", 1)[0]
        if not target or "://" in target or target.startswith("mailto:"):
            continue
        resolved = (manual.parent / target).resolve()
        assert resolved.exists(), f"broken link in {manual}: {raw_target}"


@pytest.mark.parametrize("manual", MANUALS)
def test_manual_has_no_internal_drafting_language(manual):
    text = manual.read_text(encoding="utf-8").lower()
    forbidden = ("plan §", "per plan", "phase 1", "case_00", "面向 llm")
    assert not [phrase for phrase in forbidden if phrase in text]


@pytest.mark.parametrize("manual", [ROOT / "README.md", ROOT / "README.zh-CN.md"])
def test_main_manual_is_navigable_and_covers_every_command(manual):
    text = manual.read_text(encoding="utf-8")
    assert len(text.splitlines()) < 600
    for command in (
        "mark", "skeleton", "info", "fetch", "enumerate", "shry",
        "ewald", "build",
    ):
        assert f"`{command}`" in text


def test_shry_guides_cover_each_review_stage():
    for manual in MANUALS[2:]:
        text = manual.read_text(encoding="utf-8")
        for stage in ("prepare", "count", "enum", "verify", "postprocess"):
            assert stage in text
