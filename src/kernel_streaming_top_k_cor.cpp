#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/strings.hpp>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
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

struct CorEntry {
    double cor;
    int j;
};
struct PosCmp {
    bool operator()(const CorEntry& a, const CorEntry& b) const { return a.cor > b.cor; }
};
struct NegCmp {
    bool operator()(const CorEntry& a, const CorEntry& b) const { return a.cor < b.cor; }
};

[[cpp11::register]]
data_frame streaming_top_k_cor_cpp(
    doubles_matrix<> X,
    int k,
    bool diag,
    double min_cor,
    strings row_names,
    int block_rows
) {
    const int G = (int) X.nrow();
    const int M = (int) X.ncol();
    if (k <= 0) stop("k must be positive");
    if (M < 2) stop("need at least 2 columns to compute correlations");
    if ((int) row_names.size() != G) stop("row_names length mismatch");
    if (block_rows <= 0) block_rows = 512;

    // Step 1: build a row-major copy of X, centered per row and L2-normalized
    // per row, so dot(Xn[i,:], Xn[j,:]) == Pearson correlation(X[i,:], X[j,:]).
    // Rows with <2 valid observations or zero variance yield zero rows in Xn
    // (their correlations evaluate to 0, but we mark them for NA cor below).
    std::vector<double> Xn((size_t) G * (size_t) M);
    std::vector<char> row_valid(G, 1);
    for (int i = 0; i < G; ++i) {
        double sum = 0.0;
        int n_valid = 0;
        for (int m = 0; m < M; ++m) {
            double v = X(i, m);
            if (!ISNA(v) && !std::isnan(v)) { sum += v; n_valid++; }
        }
        if (n_valid < 2) {
            row_valid[i] = 0;
            for (int m = 0; m < M; ++m) Xn[(size_t) i * M + m] = 0.0;
            continue;
        }
        double mean = sum / n_valid;
        double ss = 0.0;
        for (int m = 0; m < M; ++m) {
            double v = X(i, m);
            if (ISNA(v) || std::isnan(v)) v = mean;
            double c = v - mean;
            Xn[(size_t) i * M + m] = c;
            ss += c * c;
        }
        if (ss == 0.0) {
            row_valid[i] = 0;
            for (int m = 0; m < M; ++m) Xn[(size_t) i * M + m] = 0.0;
            continue;
        }
        double inv = 1.0 / std::sqrt(ss);
        for (int m = 0; m < M; ++m) Xn[(size_t) i * M + m] *= inv;
    }

    // Step 2: tile rows into blocks. For each block I of size up to block_rows,
    // compute C_block(bsize x G) = Xn[I,:](bsize x M) %*% t(Xn)(M x G).
    //
    // BLAS dgemm expects COLUMN-major. We trick it:
    //   dgemm("T","N", G, bsize, M, 1, Xn_full, M, Xn_block, M, 0, C_out, G)
    // computes C_out as a (G x bsize) column-major matrix. Memory-wise this is
    // identical to a (bsize x G) row-major matrix — which is exactly the view
    // we want: C_out[r * G + j] == dot(Xn[I+r, :], Xn[j, :]).

    std::vector<std::priority_queue<CorEntry, std::vector<CorEntry>, PosCmp>> pos_heaps(G);
    std::vector<std::priority_queue<CorEntry, std::vector<CorEntry>, NegCmp>> neg_heaps(G);

    std::vector<std::string> row_name_vec(G);
    for (int i = 0; i < G; ++i) row_name_vec[i] = row_names[i];

    std::vector<double> C_block;
    for (int i0 = 0; i0 < G; i0 += block_rows) {
        int i1 = std::min(i0 + block_rows, G);
        int bsize = i1 - i0;
        C_block.assign((size_t) bsize * (size_t) G, 0.0);

        {
            char trans_a = 'T', trans_b = 'N';
            int mm = G, nn = bsize, kk = M;
            double alpha = 1.0, beta = 0.0;
            int lda = M, ldb = M, ldc = G;
            F77_CALL(dgemm)(&trans_a, &trans_b, &mm, &nn, &kk, &alpha,
                            Xn.data(), &lda,
                            Xn.data() + (size_t) i0 * M, &ldb,
                            &beta, C_block.data(), &ldc FCONE FCONE);
        }

        // Update heaps for each block row in parallel.
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int r = 0; r < bsize; ++r) {
            int i = i0 + r;
            if (!row_valid[i]) continue;
            auto& ph = pos_heaps[i];
            auto& nh = neg_heaps[i];
            const double* row = C_block.data() + (size_t) r * G;
            for (int j = 0; j < G; ++j) {
                if (!row_valid[j]) continue;
                if (!diag && j == i) continue;
                double c = row[j];
                if (ISNA(c) || std::isnan(c)) continue;
                if (std::fabs(c) < min_cor) continue;
                if (c > 0) {
                    if ((int) ph.size() < k) ph.push({c, j});
                    else if (c > ph.top().cor) { ph.pop(); ph.push({c, j}); }
                } else if (c < 0) {
                    if ((int) nh.size() < k) nh.push({c, j});
                    else if (c < nh.top().cor) { nh.pop(); nh.push({c, j}); }
                }
            }
        }
    }

    // Step 3: drain heaps into long-form output.
    std::vector<std::string> out_g1, out_g2, out_type;
    std::vector<double> out_cor;
    size_t cap = (size_t) G * (size_t) 2 * (size_t) k;
    out_g1.reserve(cap); out_g2.reserve(cap);
    out_type.reserve(cap); out_cor.reserve(cap);

    for (int i = 0; i < G; ++i) {
        std::vector<CorEntry> tmp;
        while (!pos_heaps[i].empty()) { tmp.push_back(pos_heaps[i].top()); pos_heaps[i].pop(); }
        std::sort(tmp.begin(), tmp.end(),
                  [](const CorEntry& a, const CorEntry& b){ return a.cor > b.cor; });
        for (auto& e : tmp) {
            out_g1.push_back(row_name_vec[i]);
            out_g2.push_back(row_name_vec[e.j]);
            out_cor.push_back(e.cor);
            out_type.push_back("pos");
        }
        tmp.clear();
        while (!neg_heaps[i].empty()) { tmp.push_back(neg_heaps[i].top()); neg_heaps[i].pop(); }
        std::sort(tmp.begin(), tmp.end(),
                  [](const CorEntry& a, const CorEntry& b){ return a.cor < b.cor; });
        for (auto& e : tmp) {
            out_g1.push_back(row_name_vec[i]);
            out_g2.push_back(row_name_vec[e.j]);
            out_cor.push_back(e.cor);
            out_type.push_back("neg");
        }
    }

    R_xlen_t n = (R_xlen_t) out_g1.size();
    writable::strings r_g1(n), r_g2(n), r_type(n);
    writable::doubles r_cor(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        r_g1[i] = out_g1[i];
        r_g2[i] = out_g2[i];
        r_cor[i] = out_cor[i];
        r_type[i] = out_type[i];
    }
    return writable::data_frame({
        "gene1"_nm = r_g1,
        "gene2"_nm = r_g2,
        "cor"_nm = r_cor,
        "type"_nm = r_type
    });
}
