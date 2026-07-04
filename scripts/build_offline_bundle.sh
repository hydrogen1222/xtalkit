#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/build_offline_bundle.sh [--with-enumerate] [--out DIR]

Build an offline bundle for xtalkit on a networked machine.

Options:
  --with-enumerate   Include pymatgen and its dependencies for enumerate/ewald
  --out DIR          Output directory (default: /tmp/xtalkit-offline)
  -h, --help         Show this help
EOF
}

WITH_ENUMERATE=0
OUTDIR="/tmp/xtalkit-offline"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-enumerate)
      WITH_ENUMERATE=1
      shift
      ;;
    --out)
      OUTDIR="${2:-}"
      if [[ -z "$OUTDIR" ]]; then
        echo "[ERR] --out requires a directory" >&2
        exit 1
      fi
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      OUTDIR="$1"
      shift
      ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(cd "$ROOT" && python -c 'from xtalkit import __version__; print(__version__)')"
BUNDLE_DIR="$OUTDIR/xtalkit-offline-$VERSION"
TARBALL="$OUTDIR/xtalkit-offline-$VERSION.tar.gz"

cd "$ROOT"

rm -rf "$BUNDLE_DIR"
mkdir -p "$BUNDLE_DIR/dist" "$BUNDLE_DIR/wheelhouse"

python -m build --wheel --outdir "$BUNDLE_DIR/dist"

deps=(gemmi rich)
if [[ "$WITH_ENUMERATE" -eq 1 ]]; then
  deps+=(pymatgen)
fi

python -m pip download \
  --dest "$BUNDLE_DIR/wheelhouse" \
  --only-binary=:all: \
  "${deps[@]}"

cat > "$BUNDLE_DIR/install.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

python -m pip install --no-index --find-links ./wheelhouse ./dist/*.whl
EOF
chmod +x "$BUNDLE_DIR/install.sh"

cat > "$BUNDLE_DIR/README.txt" <<EOF
xtalkit offline bundle $VERSION

On the target cluster:
  1. conda activate your_environment
  2. bash install.sh
  3. xtalkit --version

If you already installed pymatgen in the conda env, this bundle still works.
If you only need the core commands, you can build the bundle without --with-enumerate.
EOF

tar -C "$OUTDIR" -czf "$TARBALL" "xtalkit-offline-$VERSION"

echo "[OK] Wrote bundle directory: $BUNDLE_DIR"
echo "[OK] Wrote tarball: $TARBALL"
