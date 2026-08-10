"""Tests for TUI module — unit tests for helper functions + smoke test."""

from unittest.mock import patch

from xtalkit.tui import (
    select_wyckoff_positions,
    select_output_formats,
    resolve_cif_path,
    _prompt_float_or_none,
    _prompt_int_or_none,
    _show_main_menu,
    run_tui,
)


def test_resolve_cif_path_exists(tmp_path):
    """relative path that exists should resolve."""
    f = tmp_path / "test.cif"
    f.write_text("data_test")
    result = resolve_cif_path(str(f))
    assert result is not None
    assert result.endswith("test.cif")


def test_resolve_cif_path_not_found():
    """Non-existent file returns None."""
    result = resolve_cif_path("./nonexistent_12345.cif")
    assert result is None


def test_resolve_cif_path_current_dir(tmp_path, monkeypatch):
    """File in current directory found by relative name."""
    f = tmp_path / "myfile.cif"
    f.write_text("data")
    monkeypatch.chdir(tmp_path)
    result = resolve_cif_path("myfile.cif")
    assert result is not None


def test_select_wyckoff_positions_all():
    """'all' should return all available letters."""
    available = ["4a", "4b", "4c", "16e", "24f", "24g", "48h"]
    with patch("builtins.input", return_value="all"):
        result = select_wyckoff_positions(available)
        assert result == available


def test_select_wyckoff_positions_subset():
    """'4a,24f' should return ['4a', '24f']."""
    with patch("builtins.input", return_value="4a,24f"):
        result = select_wyckoff_positions(["4a", "4b", "24f", "48h"])
        assert result == ["4a", "24f"]


def test_select_wyckoff_positions_invalid_retry():
    """Invalid input should prompt retry, then accept valid."""
    with patch("builtins.input", side_effect=["99z", "4a"]):
        result = select_wyckoff_positions(["4a", "4b", "4c"])
        assert result == ["4a"]


def test_select_output_formats():
    with patch("builtins.input", return_value="1"):
        result = select_output_formats()
        assert result == ["cif"]

    with patch("builtins.input", return_value="3"):
        result = select_output_formats()
        assert result == ["cif", "xyz"]


def test_show_main_menu_mentions_ewald(capsys):
    _show_main_menu()
    output = capsys.readouterr().out
    assert "Ewald" in output
    assert "Fetch DB" not in output
    assert "xtalkit" in output and "fetch" in output
    assert "SHRY" in output


def test_prompt_integer_retries_instead_of_using_dangerous_default():
    with patch("builtins.input", side_effect=["not-a-number", "0", "5"]):
        assert _prompt_int_or_none("Count", minimum=1) == 5


def test_prompt_float_rejects_nonfinite_and_out_of_range():
    with patch("builtins.input", side_effect=["nan", "0", "1.2", "0.5"]):
        assert _prompt_float_or_none(
            "Occupancy", greater_than=0.0, maximum=1.0,
        ) == 0.5


def test_prompt_numeric_blank_uses_explicit_default():
    with patch("builtins.input", return_value=""):
        assert _prompt_int_or_none("Jobs", default=1, minimum=0) == 1
    with patch("builtins.input", return_value=""):
        assert _prompt_float_or_none("Timeout") is None


def test_tui_handles_closed_input_without_traceback():
    with patch("xtalkit.tui._show_main_menu"), patch(
        "xtalkit.tui._prompt", side_effect=EOFError,
    ):
        assert run_tui() == 0


def test_tui_handles_keyboard_interrupt_without_traceback():
    with patch("xtalkit.tui._show_main_menu"), patch(
        "xtalkit.tui._prompt", side_effect=KeyboardInterrupt,
    ):
        assert run_tui() == 130
