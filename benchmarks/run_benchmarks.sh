#!/usr/bin/env bash
# ==============================================================================
# MCView DAF Benchmark Runner
# ==============================================================================
#
# Sets up the conda environment, Julia env vars, and runs the benchmark suite.
#
# Usage:
#   ./benchmarks/run_benchmarks.sh                              # defaults
#   ./benchmarks/run_benchmarks.sh --daf /path/to/daf           # custom DAF
#   ./benchmarks/run_benchmarks.sh --iterations 10              # more iterations
#   ./benchmarks/run_benchmarks.sh --benchmarks "single_gene_egc,full_egc_matrix"
#
# Environment:
#   MCVIEW_DAF_PATH   - Override default DAF path
#   MCVIEW_CONDA_ENV  - Override conda environment name (default: dafr-mcview)
#   MCVIEW_PKG_ROOT   - Override package root directory
# ==============================================================================

set -euo pipefail

# -- Defaults ------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="${MCVIEW_PKG_ROOT:-$(dirname "$SCRIPT_DIR")}"
DAF_PATH="${MCVIEW_DAF_PATH:-/home/obk/data/mcview/metacells_clean}"
CONDA_ENV="${MCVIEW_CONDA_ENV:-dafr-mcview}"
ITERATIONS=5
BENCHMARKS=""
OUTPUT=""

# -- Parse arguments -----------------------------------------------------------

PASSTHROUGH_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --daf)
            DAF_PATH="$2"
            shift 2
            ;;
        --iterations)
            ITERATIONS="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --benchmarks)
            BENCHMARKS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --daf <path>         Path to DAF directory or H5 file"
            echo "                       (default: /home/obk/data/mcview/metacells_clean)"
            echo "  --iterations <n>     Number of timed iterations (default: 5)"
            echo "  --output <path>      Output JSON file (default: auto-generated)"
            echo "  --benchmarks <list>  Comma-separated benchmark names (default: all)"
            echo "  --help               Show this help"
            echo ""
            echo "Environment variables:"
            echo "  MCVIEW_DAF_PATH      Override default DAF path"
            echo "  MCVIEW_CONDA_ENV     Override conda env name (default: dafr-mcview)"
            echo "  MCVIEW_PKG_ROOT      Override package root directory"
            exit 0
            ;;
        *)
            PASSTHROUGH_ARGS+=("$1")
            shift
            ;;
    esac
done

# -- Activate conda environment ------------------------------------------------

echo "Setting up conda environment: ${CONDA_ENV}"

# Try to initialize conda
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

# -- Build Rscript command -----------------------------------------------------

R_ARGS=(
    "--daf" "${DAF_PATH}"
    "--iterations" "${ITERATIONS}"
    "--pkg-root" "${PKG_ROOT}"
)

if [[ -n "${OUTPUT}" ]]; then
    R_ARGS+=("--output" "${OUTPUT}")
fi

if [[ -n "${BENCHMARKS}" ]]; then
    R_ARGS+=("--benchmarks" "${BENCHMARKS}")
fi

# Add any passthrough args
R_ARGS+=("${PASSTHROUGH_ARGS[@]+"${PASSTHROUGH_ARGS[@]}"}")

# -- Verify prerequisites -----------------------------------------------------

echo ""
echo "Configuration:"
echo "  Package root:  ${PKG_ROOT}"
echo "  DAF path:      ${DAF_PATH}"
echo "  Iterations:    ${ITERATIONS}"
echo "  Conda env:     ${CONDA_ENV}"
if [[ -n "${BENCHMARKS}" ]]; then
    echo "  Benchmarks:    ${BENCHMARKS}"
fi
echo ""

if [[ ! -f "${PKG_ROOT}/DESCRIPTION" ]]; then
    echo "ERROR: DESCRIPTION not found at ${PKG_ROOT}. Is --pkg-root correct?"
    exit 1
fi

if [[ ! -d "${DAF_PATH}" ]] && [[ ! -f "${DAF_PATH}" ]]; then
    echo "ERROR: DAF path not found: ${DAF_PATH}"
    exit 1
fi

# -- Run the benchmark ---------------------------------------------------------

echo "Starting benchmark..."
echo ""

Rscript "${PKG_ROOT}/benchmarks/benchmark_daf.R" "${R_ARGS[@]}"

echo ""
echo "Benchmark complete."
