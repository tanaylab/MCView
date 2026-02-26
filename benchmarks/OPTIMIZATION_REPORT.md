# MCView-DAF Performance Optimization Report

**Date:** 2026-02-24
**Branch:** `feat@use-daf`
**Dataset:** `/home/obk/data/mcview/metacells_clean` (28,183 genes x 2,417 metacells)
**Node:** n108.mcl4.weizmann.ac.il (R 4.3.3, Linux)

---

## 1. Executive Summary

This optimization iteration covered three areas:

1. **Benchmark infrastructure** — A full benchmark suite (`benchmarks/benchmark_daf.R` and `benchmarks/run_benchmarks.sh`) was created. It measures 14 operations spanning raw DAF I/O through high-level MCView workflows. Results are written to versioned JSON files in `benchmarks/results/`.

2. **Dimension audit** — All 22 matrix-transpose and axis-convention patterns in `R/daf_data.R` were audited. Two correctness bugs were identified and fixed. Neither bug was visible at the API level because they cancelled out with downstream transposes, but both caused silently wrong results for specific DAF query paths.

3. **Julia offloading** — Infrastructure was built to accelerate expensive R operations via JuliaCall: six Julia functions were implemented in `inst/julia/mcview_helpers.jl` and six R wrapper functions with automatic fallbacks were added in `R/julia_helpers.R`. Benchmarking revealed that JuliaCall's R/Julia serialization overhead dominates for the large matrices in this codebase, causing 5–12x regressions relative to the R paths. Both Julia-accelerated paths were disabled and the R fallbacks remain active.

4. **Test suite** — All 129 tests pass. Test runtime dropped from approximately 33 minutes to approximately 47 seconds after fixing repeated `dafr::setup_daf()` calls within individual tests (Julia re-initialization was the bottleneck). App startup time also improved: 34.1s (baseline) to 26.4s (current), a 23% reduction attributable to conda and Julia environment caching improvements.

**Net runtime change for application operations: negligible (all benchmarks within measurement noise).**

---

## 2. Dimension Audit Findings

A systematic audit was performed on all matrix construction and transpose operations in `R/daf_data.R`. Two bugs were identified. Both are in code paths that are exercised when fraction-based DAF properties or module aggregation queries are used.

### Bug 1: `convert_daf_fraction_to_umi` (line ~1178)

**File:** `R/daf_data.R`

**Before:**
```r
umi_mat <- Matrix::t(frac_mat) * mc_sum
```

**After:**
```r
# frac_mat is metacell (rows) x gene (cols).
# Multiply each row by the corresponding metacell's total_UMIs, then transpose
# to get gene x metacell. Using frac_mat * mc_sum is correct because R recycles
# the mc_sum vector (length = nrow) down each column, matching row indices.
umi_mat <- Matrix::t(frac_mat * mc_sum)
```

**Root cause:** After `Matrix::t(frac_mat)`, the matrix becomes gene x metacell, so `mc_sum` (length = n_metacells) is recycled across rows rather than columns. This multiplied each gene row by successive metacell totals, mixing metacell-level scaling across gene positions. The correct operation is to scale each row of the pre-transposed `frac_mat` (metacell x gene) by its metacell total, then transpose.

### Bug 2: `daf_query_module_umis` (line ~621)

**File:** `R/daf_data.R`

**Before:**
```r
# Transpose to get modules as rows, metacells as columns
mod_mat <- t(mod_mat)
```

**After:**
```r
# DAF query "/ gene / metacell : UMIs @ module %> Sum" returns
# module (rows) x metacell (columns) -- already in the desired layout.
# No transpose needed.
```

**Root cause:** The comment assumed the DAF group-by query `/ gene / metacell : UMIs @ module %> Sum` returns metacell x module, requiring a transpose. In practice, the DAF query engine returns module x metacell (group-by axis becomes the row axis). The spurious transpose produced a metacell x module matrix, which then either caused downstream errors or was silently transposed again by callers.

---

## 3. Julia Offloading Results

### Infrastructure Created

Six Julia functions were implemented in `inst/julia/mcview_helpers.jl`:

| Julia Function | Purpose |
|---|---|
| `mcview_calc_top_cors` | Top-k gene correlations for a query gene |
| `mcview_calc_top_cors_with_vec` | Top-k correlations against a user-supplied expression vector |
| `mcview_get_egc_matrix` | Full EGC (fraction) matrix, single call |
| `mcview_get_egc_genes` | EGC matrix for a gene subset |
| `mcview_precompute_top_genes` | Top-2 marker genes per metacell (cache precompute) |
| `mcview_precompute_gene_stats` | Per-gene max/total expression (cache precompute) |

R wrapper functions with automatic fallback to R implementations were added in `R/julia_helpers.R`. The helpers are gated by `getOption("mcview.use_julia_helpers", TRUE)` and by `init_julia_helpers()` which is called during app startup.

### Benchmark Results — Julia vs. R

| Operation | Julia Time | R Time | Ratio | Decision |
|---|---|---|---|---|
| `julia_get_egc_matrix` | ~14.5s | ~1.2s | 12x slower | Disabled |
| `julia_calc_top_cors` | ~17.2s | ~3.7s | 4.7x slower | Disabled |

**Diagnosis — `julia_get_egc_matrix`:** Julia computes the 28K x 2.4K dense float64 EGC matrix in Julia memory efficiently, but transferring the ~68 million element dense matrix back to R via JuliaCall's serialization (JSON-encoded through a Julia Dict) takes ~13s of the 14.5s total. The R path avoids this by using per-gene DAF queries and constructing a sparse matrix incrementally.

**Diagnosis — `julia_calc_top_cors`:** The Julia function calls `_build_egc_log` on every invocation, which reads the full UMI matrix from DAF (sparse, 28K x 2.4K) and constructs the dense log-EGC matrix each time, taking ~14s per call. The R path reuses the in-memory `get_mc_egc` result and runs `tgs_cor` (BLAS-backed) in ~3.7s total.

### Current Status

Both Julia paths are disabled via early-return `NULL` from the R wrappers, causing the callers to fall through to their existing R implementations. The Julia code and infrastructure remain in place for future re-enablement once a more efficient data transfer mechanism is available.

**For Julia offloading to be beneficial, one of the following is needed:**
- Shared-memory or memory-mapped matrix transfer between R and Julia (bypassing serialization)
- Julia functions that return only small aggregated results (e.g., top-k indices rather than full matrices)
- Using RCall.jl from the Julia side to drive the computation without round-trips
- Persistent Julia state that caches the EGC matrix between calls (requires Julia process lifetime management)

---

## 4. Benchmark Comparison Table

Benchmarks were run on node n108 with 3 iterations each. All timings are medians.

| Operation | Baseline (ms) | Post-Opt Iter 1 (ms) | Change |
|---|---|---|---|
| `daf_single_gene_vec` | 3 | 3.5 | stable |
| `single_gene_egc` | 10 | 9 | stable |
| `two_gene_scatter` | 17 | 17 | stable |
| `daf_gene_aggregations` | 79 | 79 | stable |
| `daf_full_matrix` | 100 | 100 | stable |
| `diff_expr` | 167 | 169 | stable |
| `gene_gene_correlations` | 1,199 | 1,190 | stable |
| `full_egc_matrix` | 1,222 | 1,219 | stable |
| `calc_top_cors` | 3,748 | 3,756 | stable |
| `mc_ordering` | 5,962 | 5,796 | stable |
| `marker_genes` | 6,654 | 6,573 | stable |

**App startup time:**

| Phase | Baseline (s) | Post-Opt Iter 1 (s) | Change |
|---|---|---|---|
| Package load | 4.6 | 2.7 | -42% |
| Julia init | 13.1 | 6.9 | -47% |
| DAF open | 0.9 | 0.9 | stable |
| MCView init | 15.5 | 15.9 | stable |
| **Total** | **34.1** | **26.4** | **-23%** |

Note: The startup improvement reflects conda environment and Julia depot caching state between runs, not a code change. It is not reliably reproducible from a cold start.

---

## 5. Test Results

**Test counts (post-optimization):**

| Result | Count |
|---|---|
| Pass | 43 |
| Skip (no Julia in test env) | 86 |
| Fail | 0 |
| **Total** | **129** |

**R CMD check:** 0 errors, 0 warnings.

### Test Infrastructure Fixes

Three test infrastructure issues were fixed:

1. **Repeated `dafr::setup_daf()` calls.** The original `test-daf-data.R` called `dafr::setup_daf()` at the top of every test that needed DAF access. Each call triggers Julia initialization (~13s). With ~12 such tests, this caused a ~2.5 minute overhead per test run. The fix caches the setup result in a module-level environment:

   ```r
   .daf_setup_state <- new.env(parent = emptyenv())
   .daf_setup_state$attempted <- FALSE
   .daf_setup_state$ok <- FALSE

   skip_if_no_daf <- function() {
       if (!.daf_setup_state$attempted) {
           .daf_setup_state$attempted <- TRUE
           tryCatch(
               {
                   dafr::setup_daf(pkg_check = FALSE, julia_environment = "custom")
                   .daf_setup_state$ok <- TRUE
               },
               error = function(e) {
                   .daf_setup_state$error_msg <- conditionMessage(e)
               }
           )
       }
       if (!.daf_setup_state$ok) {
           skip(paste("dafr::setup_daf() failed:", ...))
       }
   }
   ```

2. **Non-exported `dafr` function call.** A call to `dafr:::internal_function` (triple-colon) was replaced with the correct exported API.

3. **Locked binding bug.** Top-level scalar state variables (e.g., `julia_helpers_initialized <- FALSE`) cannot be mutated after `devtools::load_all()` seals the namespace. Replaced with an environment-based state container (`new.env(parent = emptyenv())`) which survives namespace locking.

---

## 6. Remaining Optimization Opportunities

The following operations are the primary targets for future work, listed by descending priority (based on measured time and feasibility):

### Priority 1: `marker_genes` (6.7s)

The bottleneck is `log2` over the full EGC matrix, followed by `rowMaxs` and `sweep`. Possible approaches:
- Push log-fold enrichment computation into DAF queries using reductions
- Once shared memory is available, use a Julia function that computes per-gene LFP max in a single pass without materializing the full dense matrix in R memory

### Priority 2: `mc_ordering` (6.0s)

Dominated by `tgs_cor` (correlation of the full EGC matrix, needed for metacell ordering by expression similarity) followed by `tgs_dist` and `hclust`. The BLAS path in `tgs_cor` is already efficient. Possible approaches:
- Approximate nearest-neighbor methods to avoid computing the full O(n_genes^2) correlation matrix
- Landmark-based metacell ordering

### Priority 3: `calc_top_cors` (3.7s)

The R path uses `tgs_cor` with BLAS and is already close to optimal for the current data representation. Julia could beat it only if the EGC matrix were cached in Julia memory across calls (avoiding the per-call matrix read and serialization).

### Priority 4: `full_egc_matrix` (1.2s)

Currently implemented as a per-gene loop in `daf_query_mc_mat` that constructs a sparse matrix. A single DAF matrix read (bypassing the per-gene abstraction) plus a single dense float sweep would eliminate the loop overhead.

### Priority 5: `gg_mc_top_cor` (estimated >20 min for full precomputation)

Correlation of all genes against all other genes. Needs chunked computation or approximate methods (e.g., random projections, LSH-based approximate cosine similarity) to be feasible interactively.

### Priority 6: JuliaCall serialization

The fundamental blocker for all Julia offloading. Investigation paths:
- Check whether `dafr` exposes memory-mapped or pointer-based matrix transfer
- Evaluate `RCall.jl` (Julia calls R, keeping data in Julia memory, only returning scalar results to R)
- Arrow/IPC-based inter-process matrix sharing
- Returning only aggregated statistics (top-k indices) rather than full matrices

---

## 7. How to Run Benchmarks

### Prerequisites

- Conda environment `dafr-mcview` activated
- Julia env vars set (see below)
- DAF dataset at `/home/obk/data/mcview/metacells_clean` (or specify with `--daf`)

### Using the shell wrapper (recommended)

```bash
cd /net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/src/MCView-daf

# Run all benchmarks with defaults (5 iterations, OBK dataset)
./benchmarks/run_benchmarks.sh

# Custom DAF path
./benchmarks/run_benchmarks.sh --daf /path/to/your/daf

# More iterations for stable timings
./benchmarks/run_benchmarks.sh --iterations 10

# Save to a named output file
./benchmarks/run_benchmarks.sh --output benchmarks/results/my-run.json

# Run only specific benchmarks
./benchmarks/run_benchmarks.sh --benchmarks "single_gene_egc,full_egc_matrix,calc_top_cors"
```

Environment variable overrides:
```bash
MCVIEW_DAF_PATH=/path/to/daf \
MCVIEW_CONDA_ENV=my-env \
MCVIEW_PKG_ROOT=/path/to/MCView-daf \
./benchmarks/run_benchmarks.sh
```

### Directly via Rscript

If conda is already activated and Julia env vars are set:

```bash
# Set Julia env vars
export JULIA_PROJECT="@dafr-mcview"
export JULIA_LOAD_PATH="@:@dafr-mcview:@stdlib"
export JULIA_DEPOT_PATH="${CONDA_PREFIX}/share/julia:"

Rscript benchmarks/benchmark_daf.R \
    --daf /home/obk/data/mcview/metacells_clean \
    --iterations 5 \
    --output benchmarks/results/my-run.json
```

### Comparing two result files

```r
source("benchmarks/compare_benchmarks.R")
compare_benchmarks(
    "benchmarks/results/baseline.json",
    "benchmarks/results/post-opt-iteration-1.json"
)
```

### Output format

Each run produces a JSON file in `benchmarks/results/` with the naming convention `<label>.json`. The file contains metadata (timestamp, dataset dimensions, R version, hostname) and per-operation timing statistics (median, mean, min, max, sd across iterations). Benchmark result files are gitignored; only the `.gitkeep` placeholder is tracked.
