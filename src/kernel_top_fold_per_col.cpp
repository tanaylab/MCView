#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/strings.hpp>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace cpp11;
using namespace cpp11::literals;

// Per-column bounded top-K of fold-point values with MCView marker thresholds.
// Applied in order, matching select_top_fold_genes_per_metacell:
//   1. add fold_reg to every non-NA value
//   2. mask values < min_val to NA (skip)
//   3. rank by abs(x) when use_abs else raw x, descending
//   4. return k rows per column; rows without enough qualifying genes get
//      NA gene / NA fp but keep metacell + rank.
//
// Returns data.frame(metacell, gene, rank, fp) length = k * ncol(fp).

struct HeapEntry {
    double score;  // comparison key (abs(x) or x)
    double value;  // signed value to return
    int row;       // gene index
};

// Min-heap by score: heap.top() is the smallest score in the heap, evicted
// when a larger score is offered.
struct HeapCmp {
    bool operator()(const HeapEntry& a, const HeapEntry& b) const {
        return a.score > b.score;
    }
};

[[cpp11::register]]
data_frame top_fold_per_col_cpp(
    doubles_matrix<> fp,
    int k,
    double min_val,
    double fold_reg,
    bool use_abs,
    strings row_names,
    strings col_names
) {
    const R_xlen_t n_rows = fp.nrow();
    const R_xlen_t n_cols = fp.ncol();

    if (k <= 0) stop("k must be positive");
    if ((R_xlen_t) row_names.size() != n_rows)
        stop("row_names length does not match fp rows");
    if ((R_xlen_t) col_names.size() != n_cols)
        stop("col_names length does not match fp cols");

    std::vector<std::string> row_name_vec(n_rows);
    for (R_xlen_t i = 0; i < n_rows; ++i) row_name_vec[i] = row_names[i];
    std::vector<std::string> col_name_vec(n_cols);
    for (R_xlen_t j = 0; j < n_cols; ++j) col_name_vec[j] = col_names[j];

    R_xlen_t n_out = (R_xlen_t) k * n_cols;
    std::vector<std::string> out_metacell(n_out);
    std::vector<std::string> out_gene(n_out);
    std::vector<int> out_rank(n_out);
    std::vector<double> out_fp(n_out);
    std::vector<char> gene_is_na(n_out, 0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (R_xlen_t j = 0; j < n_cols; ++j) {
        std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapCmp> heap;

        for (R_xlen_t i = 0; i < n_rows; ++i) {
            double v = fp(i, j);
            if (ISNA(v) || std::isnan(v)) continue;
            v += fold_reg;
            if (v < min_val) continue;
            double score = use_abs ? std::fabs(v) : v;

            if ((int) heap.size() < k) {
                heap.push({score, v, (int) i});
            } else if (score > heap.top().score) {
                heap.pop();
                heap.push({score, v, (int) i});
            }
        }

        std::vector<HeapEntry> ordered;
        ordered.reserve(heap.size());
        while (!heap.empty()) {
            ordered.push_back(heap.top());
            heap.pop();
        }
        std::reverse(ordered.begin(), ordered.end());

        for (int r = 0; r < k; ++r) {
            R_xlen_t out_idx = j * k + r;
            out_metacell[out_idx] = col_name_vec[j];
            out_rank[out_idx] = r + 1;
            if (r < (int) ordered.size()) {
                out_gene[out_idx] = row_name_vec[ordered[r].row];
                out_fp[out_idx] = ordered[r].value;
            } else {
                out_gene[out_idx] = "";
                gene_is_na[out_idx] = 1;
                out_fp[out_idx] = NA_REAL;
            }
        }
    }

    writable::strings r_metacell(n_out);
    writable::strings r_gene(n_out);
    writable::integers r_rank(n_out);
    writable::doubles r_fp(n_out);
    for (R_xlen_t i = 0; i < n_out; ++i) {
        r_metacell[i] = out_metacell[i];
        if (gene_is_na[i]) {
            r_gene[i] = NA_STRING;
        } else {
            r_gene[i] = out_gene[i];
        }
        r_rank[i] = out_rank[i];
        r_fp[i] = out_fp[i];
    }

    return writable::data_frame({
        "metacell"_nm = r_metacell,
        "gene"_nm = r_gene,
        "rank"_nm = r_rank,
        "fp"_nm = r_fp
    });
}
