#!/usr/bin/env bash
#
# build_enumlib.sh — compile Hart & Forcade's enumlib (msg-byu/enumlib) from
# source with the system Fortran compiler, then install the `enum.x` and
# `makestr.x` executables into a user-local directory.
#
# No root required. The binaries are picked up automatically by
# `xtalkit._env.setup_for_enumlib()` (via PATH), so after running this once,
# `xtalkit enumerate` works without any manual PATH configuration.
#
# Overrides:
#   F90                  Fortran compiler (default: gfortran)
#   XTALKIT_ENUMLIB_BIN  install directory (default: ~/.local/share/xtalkit/bin)
#
set -euo pipefail

F90="${F90:-gfortran}"
PREFIX="${XTALKIT_ENUMLIB_BIN:-$HOME/.local/share/xtalkit/bin}"

# --- prerequisite checks -----------------------------------------------------
for cmd in "$F90" git make; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERR] '$cmd' not found in PATH." >&2
        case "$cmd" in
            gfortran|"$F90") echo "      Install a Fortran compiler, e.g.:" >&2
                            echo "        Debian/Ubuntu: sudo apt install gfortran" >&2
                            echo "        Fedora:        sudo dnf install gcc-gfortran" >&2
                            echo "        macOS:         brew install gcc" >&2 ;;
            git) echo "      Install git (e.g. sudo apt install git)." >&2 ;;
            make) echo "      Install build tools (e.g. sudo apt install make)." >&2 ;;
        esac
        exit 1
    fi
done

echo "[1/4] Compiler: $F90   |   Install prefix: $PREFIX"

# --- clone (full clone so `git describe` works for the version tag) ---------
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
echo "[2/4] Cloning enumlib (+ symlib submodule) into $WORK/enumlib ..."
git clone --recursive https://github.com/msg-byu/enumlib.git "$WORK/enumlib"

# --- build symlib (dependency) then enum.x + makestr.x ----------------------
echo "[3/4] Compiling symlib ..."
( cd "$WORK/enumlib/symlib/src" && make "F90=$F90" )

echo "[4/4] Compiling enum.x + makestr.x ..."
( cd "$WORK/enumlib/src" && make "F90=$F90" enum.x makestr.x )

# --- install -----------------------------------------------------------------
mkdir -p "$PREFIX"
cp "$WORK/enumlib/src/enum.x" "$WORK/enumlib/src/makestr.x" "$PREFIX/"
chmod 0755 "$PREFIX/enum.x" "$PREFIX/makestr.x"

cat <<EOF

[OK] Installed:
      $PREFIX/enum.x
      $PREFIX/makestr.x

xtalkit discovers these automatically — no PATH setup needed. To use a
different location, set XTALKIT_ENUMLIB_BIN and re-run this script.
EOF
