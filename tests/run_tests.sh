#!/usr/bin/env bash
# ==============================================================================
# MCView Test Runner with Julia Environment
# ==============================================================================
#
# Sets up the conda environment, Julia env vars, and runs devtools::test().
# This ensures tests use the correct Julia binary (conda's) instead of
# whatever is on PATH (e.g., juliaup's version).
#
# Usage:
#   ./tests/run_tests.sh                              # run all tests
#   ./tests/run_tests.sh --filter "daf-data"          # run specific test file
#
# Environment:
#   MCVIEW_CONDA_ENV  - Override conda environment name (default: dafr-mcview)
#   MCVIEW_PKG_ROOT   - Override package root directory
# ==============================================================================

set -euo pipefail

# -- Defaults ------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="${MCVIEW_PKG_ROOT:-$(dirname "$SCRIPT_DIR")}"
CONDA_ENV="${MCVIEW_CONDA_ENV:-dafr-mcview}"
FILTER=""

# -- Parse arguments -----------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --filter)
            FILTER="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --filter <pattern>   Filter test files by pattern (passed to devtools::test)"
            echo "  --help               Show this help"
            echo ""
            echo "Environment variables:"
            echo "  MCVIEW_CONDA_ENV     Override conda env name (default: dafr-mcview)"
            echo "  MCVIEW_PKG_ROOT      Override package root directory"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

# -- Activate conda environment ------------------------------------------------

echo "Setting up conda environment: ${CONDA_ENV}"

if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook 2>/dev/null)" || true
fi

conda activate "${CONDA_ENV}" 2>/dev/null || {
    echo "ERROR: Could not activate conda environment '${CONDA_ENV}'."
    echo "Make sure conda is initialized and the environment exists."
    echo "Run: conda env list"
    exit 1
}

echo "Conda prefix: ${CONDA_PREFIX}"

# -- Set Julia environment variables -------------------------------------------

export JULIA_PROJECT="@dafr-mcview"
export JULIA_LOAD_PATH="@:@dafr-mcview:@stdlib"
export JULIA_DEPOT_PATH="${CONDA_PREFIX}/share/julia:"

# Use sysimage if available
SYSIMAGE="${CONDA_PREFIX}/share/julia/sysimage_daf.so"
if [[ -f "${SYSIMAGE}" ]]; then
    export JULIA_SYSIMAGE="${SYSIMAGE}"
    echo "Using Julia sysimage: ${SYSIMAGE}"
fi

# -- Verify prerequisites -----------------------------------------------------

echo ""
echo "Configuration:"
echo "  Package root:  ${PKG_ROOT}"
echo "  Conda env:     ${CONDA_ENV}"
echo "  Julia binary:  $(which julia 2>/dev/null || echo 'NOT FOUND')"
echo "  Julia version: $(julia --version 2>/dev/null || echo 'N/A')"
if [[ -n "${FILTER}" ]]; then
    echo "  Filter:        ${FILTER}"
fi
echo ""

if [[ ! -f "${PKG_ROOT}/DESCRIPTION" ]]; then
    echo "ERROR: DESCRIPTION not found at ${PKG_ROOT}. Is --pkg-root correct?"
    exit 1
fi

# -- Run tests -----------------------------------------------------------------

echo "Starting tests..."
echo ""

if [[ -n "${FILTER}" ]]; then
    Rscript -e "devtools::test('${PKG_ROOT}', filter = '${FILTER}')"
else
    Rscript -e "devtools::test('${PKG_ROOT}')"
fi

echo ""
echo "Tests complete."
