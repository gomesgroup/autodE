#!/bin/bash
# Setup script for autodE 1.4.5.gpggrp.3 (gomesgroup fork with ORCA 6.x + RACE-TS)
# Auto-detects architecture (x86_64 or ARM64)

AUTODE_ROOT="/mnt/beegfs/software/autode"

# Detect architecture
ARCH=$(uname -m)
case "$ARCH" in
    x86_64)
        ENV_DIR="${AUTODE_ROOT}/env"
        ;;
    aarch64)
        ENV_DIR="${AUTODE_ROOT}/env-arm64"
        ;;
    *)
        echo "Error: Unsupported architecture: $ARCH"
        return 1
        ;;
esac

# Check if environment exists
if [[ ! -d "$ENV_DIR" ]]; then
    echo "Error: autodE environment not found at $ENV_DIR"
    return 1
fi

# Set up environment (conda-style)
export PATH="${ENV_DIR}/bin:$PATH"
export LD_LIBRARY_PATH="${ENV_DIR}/lib:$LD_LIBRARY_PATH"
export AUTODE_HOME="$AUTODE_ROOT"

# Verify autodE is importable
if ! "${ENV_DIR}/bin/python" -c "import autode" 2>/dev/null; then
    echo "Warning: autodE module not importable - check installation"
fi

echo "autodE 1.4.5.gpggrp.3 loaded ($ARCH)"
echo ""
echo "Features: ORCA 6.x integration, RACE-TS conformer generation"
echo ""
echo "Quick start:"
echo "  python -c 'import autode as ade; print(ade.__version__)'"
echo ""
echo "For ORCA 6.x, also run:"
echo "  source /mnt/beegfs/software/orca-6.1.1/setup-orca.sh"
echo "  source /mnt/beegfs/software/xtb-6.7.1/setup-xtb.sh"
