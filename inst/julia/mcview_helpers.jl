# mcview_helpers.jl - Julia helper functions for MCView performance-critical paths
#
# These functions operate directly on DAF data in Julia, avoiding expensive
# round-trips through R for large matrix operations. They are loaded via
# JuliaCall::julia_source() from R/julia_helpers.R.
#
# Requirements:
#   - DataAxesFormats (loaded by dafr)
#   - LinearAlgebra (stdlib, loaded by dafr)
#   - Statistics (stdlib)
#   - SparseArrays (stdlib)

using DataAxesFormats
using Statistics
using LinearAlgebra
using SparseArrays

# ==============================================================================
# EGC Matrix Cache
# ==============================================================================

# Cache the EGC log matrix to avoid rebuilding on every calc_top_cors call.
# The cache is keyed by the DAF object identity (pointer) and epsilon value.
const _EGC_CACHE = Ref{Any}(nothing)

struct EGCCache
    daf_id::UInt64       # objectid of the DAF reader
    egc_epsilon::Float64
    egc_log::Matrix{Float64}   # n_genes x n_mc, log2(fraction + eps)
    row_norms::Vector{Float64} # L2 norms of centered rows
    gene_names::Vector{String}
    metacell_names::Vector{String}
end

function _get_cached_egc(daf, egc_epsilon::Float64)
    cache = _EGC_CACHE[]
    daf_id = objectid(daf)
    if cache !== nothing && cache.daf_id == daf_id && cache.egc_epsilon == egc_epsilon
        return cache
    end

    # Build fresh cache
    egc_log, gene_names, metacell_names = _build_egc_log(daf, egc_epsilon)
    row_norms = _center_rows!(egc_log)

    cache = EGCCache(daf_id, egc_epsilon, egc_log, row_norms,
                     String.(gene_names), String.(metacell_names))
    _EGC_CACHE[] = cache
    return cache
end

# ==============================================================================
# Internal helpers
# ==============================================================================

"""
Build the EGC (fraction) matrix in Julia from UMIs + total_UMIs.
Returns (egc_log, gene_names, metacell_names) where egc_log is
n_genes x n_mc with log2(fraction + epsilon) values.

Uses column-efficient access: iterates genes (columns of UMI CSC matrix)
in the outer loop to avoid cache-unfriendly row striding.
"""
function _build_egc_log(daf, egc_epsilon::Float64)
    metacell_names = axis_vector(daf, "metacell")
    gene_names = axis_vector(daf, "gene")
    n_mc = length(metacell_names)
    n_genes = length(gene_names)

    # Get UMIs as mc x gene (CSC: columns = genes, fast column access)
    umis = get_matrix(daf, "metacell", "gene", "UMIs")
    total_umis = get_vector(daf, "metacell", "total_UMIs")

    # Precompute inverse totals
    inv_totals = Vector{Float64}(undef, n_mc)
    @inbounds for mc in 1:n_mc
        inv_totals[mc] = 1.0 / Float64(total_umis[mc])
    end

    # Build gene x mc dense matrix of log2(fraction + eps)
    # Iterate genes in OUTER loop: for CSC mc×gene matrix,
    # umis[:, g] (all metacells for gene g) is a contiguous column access
    egc_log = Matrix{Float64}(undef, n_genes, n_mc)
    @inbounds for g in 1:n_genes
        col = umis[:, g]  # fast CSC column extraction
        for mc in 1:n_mc
            egc_log[g, mc] = log2(Float64(col[mc]) * inv_totals[mc] + egc_epsilon)
        end
    end

    return egc_log, gene_names, metacell_names
end

"""
Build the raw EGC (fraction) matrix without log transform.
Returns (egc, gene_names, metacell_names) where egc is n_genes x n_mc.
"""
function _build_egc_raw(daf)
    metacell_names = axis_vector(daf, "metacell")
    gene_names = axis_vector(daf, "gene")
    n_mc = length(metacell_names)
    n_genes = length(gene_names)

    umis = get_matrix(daf, "metacell", "gene", "UMIs")
    total_umis = get_vector(daf, "metacell", "total_UMIs")

    inv_totals = Vector{Float64}(undef, n_mc)
    @inbounds for mc in 1:n_mc
        inv_totals[mc] = 1.0 / Float64(total_umis[mc])
    end

    # Iterate genes in outer loop for efficient CSC column access
    egc = Matrix{Float64}(undef, n_genes, n_mc)
    @inbounds for g in 1:n_genes
        col = umis[:, g]
        for mc in 1:n_mc
            egc[g, mc] = Float64(col[mc]) * inv_totals[mc]
        end
    end

    return egc, gene_names, metacell_names
end

"""
Compute Pearson correlation between a target vector and all rows of a
centered matrix. Uses BLAS gemv for the dot products.

center_mat: n_genes x n_mc (already row-centered, rows have mean 0)
row_norms:  n_genes vector of row L2 norms
target:     n_mc vector (already centered)
target_norm: L2 norm of target

Returns: n_genes vector of correlations
"""
function _cors_against_all(center_mat::Matrix{Float64},
                           row_norms::Vector{Float64},
                           target::Vector{Float64},
                           target_norm::Float64)
    n_genes = size(center_mat, 1)
    # BLAS gemv: dots = center_mat * target (n_genes dot products at once)
    dots = center_mat * target

    cors = Vector{Float64}(undef, n_genes)
    inv_target_norm = 1.0 / target_norm
    @inbounds for g in 1:n_genes
        if row_norms[g] < 1e-15
            cors[g] = 0.0
        else
            cors[g] = dots[g] / (row_norms[g] * target_norm)
        end
    end
    return cors
end

"""
Row-center a matrix in-place and return row norms.
Iterates column-first for cache efficiency on column-major Julia matrices.
"""
function _center_rows!(mat::Matrix{Float64})
    n_rows, n_cols = size(mat)

    # Phase 1: compute row means (accumulate column-by-column for cache efficiency)
    row_means = zeros(Float64, n_rows)
    @inbounds for c in 1:n_cols
        for r in 1:n_rows
            row_means[r] += mat[r, c]
        end
    end
    inv_cols = 1.0 / n_cols
    @inbounds for r in 1:n_rows
        row_means[r] *= inv_cols
    end

    # Phase 2: subtract means and accumulate norms (column-by-column)
    row_norm_sq = zeros(Float64, n_rows)
    @inbounds for c in 1:n_cols
        for r in 1:n_rows
            mat[r, c] -= row_means[r]
            row_norm_sq[r] += mat[r, c]^2
        end
    end

    row_norms = Vector{Float64}(undef, n_rows)
    @inbounds for r in 1:n_rows
        row_norms[r] = sqrt(row_norm_sq[r])
    end
    return row_norms
end

"""
Extract top-k from correlation vector, skipping excluded genes.
Returns (genes, cors, types) vectors.
"""
function _extract_topk(cors::Vector{Float64}, gene_names::Vector{String},
                       k::Int, exclude_set::Set{String})
    result_genes = String[]
    result_cors = Float64[]
    result_types = String[]

    # Top-k positive
    sorted_pos = sortperm(cors, rev=true)
    count = 0
    for idx in sorted_pos
        count >= k && break
        gname = gene_names[idx]
        if !(gname in exclude_set)
            push!(result_genes, gname)
            push!(result_cors, cors[idx])
            push!(result_types, "pos")
            count += 1
        end
    end

    # Top-k negative
    sorted_neg = sortperm(cors)
    count = 0
    for idx in sorted_neg
        count >= k && break
        gname = gene_names[idx]
        if !(gname in exclude_set)
            push!(result_genes, gname)
            push!(result_cors, cors[idx])
            push!(result_types, "neg")
            count += 1
        end
    end

    return result_genes, result_cors, result_types
end


# ==============================================================================
# 1. calc_top_cors - Single gene vs all correlation (replaces 3.7s R path)
# ==============================================================================

"""
    mcview_calc_top_cors(daf, gene::AbstractString,
                         egc_epsilon::Float64, k::Int64)

Compute the top-k positively and negatively correlated genes with `gene`,
using the cached centered EGC log matrix. Returns a Dict with keys
"genes", "cors", "types" as parallel vectors.

Strategy: Use cached log2(EGC + eps) matrix (centered, with precomputed norms),
then BLAS gemv for all dot products at once. Only returns small top-k result.
"""
function mcview_calc_top_cors(daf, gene::AbstractString,
                              egc_epsilon::Float64, k::Int64)
    cache = _get_cached_egc(daf, egc_epsilon)

    # Find gene index
    gene_idx = findfirst(==(gene), cache.gene_names)
    if gene_idx === nothing
        error("Gene '$gene' not found in DAF")
    end

    target = cache.egc_log[gene_idx, :]
    target_norm = cache.row_norms[gene_idx]

    if target_norm < 1e-15
        return Dict("genes" => String[], "cors" => Float64[], "types" => String[])
    end

    # BLAS-accelerated correlation against all genes
    cors = _cors_against_all(cache.egc_log, cache.row_norms, target, target_norm)

    # Extract top-k, excluding self
    exclude_set = Set([String(gene)])
    genes, cor_vals, types = _extract_topk(cors, cache.gene_names, k, exclude_set)

    return Dict("genes" => genes, "cors" => cor_vals, "types" => types)
end


"""
    mcview_calc_top_cors_with_vec(daf, data_vec::Vector{Float64},
                                  egc_epsilon::Float64, k::Int64,
                                  exclude::Vector{String})

Compute the top-k positively and negatively correlated genes with a
user-supplied expression vector. Returns same Dict structure.
"""
function mcview_calc_top_cors_with_vec(daf, data_vec::Vector{Float64},
                                       egc_epsilon::Float64, k::Int64,
                                       exclude::Vector{String})
    cache = _get_cached_egc(daf, egc_epsilon)
    n_mc = size(cache.egc_log, 2)

    if length(data_vec) != n_mc
        error("data_vec length ($(length(data_vec))) != number of metacells ($n_mc)")
    end

    # Prepare target vector (log2 + center)
    target_log = log2.(data_vec .+ egc_epsilon)
    target_mean = mean(target_log)
    target = target_log .- target_mean
    target_norm = sqrt(sum(target .^ 2))

    if target_norm < 1e-15
        return Dict("genes" => String[], "cors" => Float64[], "types" => String[])
    end

    cors = _cors_against_all(cache.egc_log, cache.row_norms, target, target_norm)

    exclude_set = Set(String.(exclude))
    genes, cor_vals, types = _extract_topk(cors, cache.gene_names, k, exclude_set)

    return Dict("genes" => genes, "cors" => cor_vals, "types" => types)
end


# ==============================================================================
# 2. Full EGC matrix retrieval (replaces 1.2s R path with per-gene queries)
# ==============================================================================

"""
    mcview_get_egc_matrix(daf)

Retrieve the full EGC matrix (gene x metacell) as a dense Float64 matrix.
Returns a Dict with "matrix", "gene_names", "metacell_names".
"""
function mcview_get_egc_matrix(daf)
    egc, gene_names, metacell_names = _build_egc_raw(daf)
    return Dict(
        "matrix" => egc,
        "gene_names" => String.(gene_names),
        "metacell_names" => String.(metacell_names)
    )
end


"""
    mcview_get_egc_genes(daf, genes::Vector{String})

Retrieve EGC for a subset of genes.
"""
function mcview_get_egc_genes(daf, genes::Vector{String})
    metacell_names = axis_vector(daf, "metacell")
    all_gene_names = axis_vector(daf, "gene")
    n_mc = length(metacell_names)

    umis = get_matrix(daf, "metacell", "gene", "UMIs")
    total_umis = get_vector(daf, "metacell", "total_UMIs")

    inv_totals = Vector{Float64}(undef, n_mc)
    @inbounds for mc in 1:n_mc
        inv_totals[mc] = 1.0 / Float64(total_umis[mc])
    end

    # Build gene index map
    gene_index_map = Dict(String(g) => i for (i, g) in enumerate(all_gene_names))
    valid_genes = String[]
    valid_indices = Int[]
    for g in genes
        idx = get(gene_index_map, g, nothing)
        if idx !== nothing
            push!(valid_genes, g)
            push!(valid_indices, idx)
        end
    end

    n_valid = length(valid_genes)
    egc = Matrix{Float64}(undef, n_valid, n_mc)

    # Iterate valid genes in outer loop for CSC column access
    @inbounds for (row, g_idx) in enumerate(valid_indices)
        col = umis[:, g_idx]
        for mc in 1:n_mc
            egc[row, mc] = Float64(col[mc]) * inv_totals[mc]
        end
    end

    return Dict(
        "matrix" => egc,
        "gene_names" => valid_genes,
        "metacell_names" => String.(metacell_names)
    )
end


# ==============================================================================
# 3. Precompute cache reductions (top genes, gene stats)
# ==============================================================================

"""
    mcview_precompute_top_genes(daf, egc_epsilon::Float64)

Compute top-2 genes per metacell by log-fold enrichment.
"""
function mcview_precompute_top_genes(daf, egc_epsilon::Float64)
    egc, gene_names, metacell_names = _build_egc_raw(daf)
    n_genes, n_mc = size(egc)

    # Compute gene means
    gene_means = vec(mean(egc, dims=2))
    @inbounds for g in 1:n_genes
        if gene_means[g] == 0.0
            gene_means[g] = 1e-10
        end
    end

    top1_gene = Vector{String}(undef, n_mc)
    top2_gene = Vector{String}(undef, n_mc)
    top1_lfp = Vector{Float64}(undef, n_mc)
    top2_lfp = Vector{Float64}(undef, n_mc)

    # egc is gene x mc (column-major), iterating mc in outer loop
    # accesses egc[:, mc] which is a contiguous column - FAST
    @inbounds for mc in 1:n_mc
        best_val = -Inf
        best_idx = 1
        second_val = -Inf
        second_idx = 1

        for g in 1:n_genes
            lfp_val = log2(egc[g, mc] / gene_means[g] + egc_epsilon)
            if lfp_val > best_val
                second_val = best_val
                second_idx = best_idx
                best_val = lfp_val
                best_idx = g
            elseif lfp_val > second_val
                second_val = lfp_val
                second_idx = g
            end
        end

        top1_gene[mc] = String(gene_names[best_idx])
        top1_lfp[mc] = best_val
        top2_gene[mc] = String(gene_names[second_idx])
        top2_lfp[mc] = second_val
    end

    return Dict(
        "top1_gene" => top1_gene,
        "top2_gene" => top2_gene,
        "top1_lfp" => top1_lfp,
        "top2_lfp" => top2_lfp
    )
end


"""
    mcview_precompute_gene_stats(daf)

Compute per-gene statistics: max_expr and total_expr across metacells.
"""
function mcview_precompute_gene_stats(daf)
    egc, gene_names, _ = _build_egc_raw(daf)

    max_expr = vec(maximum(egc, dims=2))
    total_expr = vec(sum(egc, dims=2))

    return Dict(
        "max_expr" => max_expr,
        "total_expr" => total_expr,
        "gene_names" => String.(gene_names)
    )
end


# ==============================================================================
# 4. calc_gg_mc_top_cor - All-vs-all gene correlation top-k (replaces 572s R path)
# ==============================================================================

"""
    mcview_calc_gg_mc_top_cor(daf, egc_epsilon::Float64, k::Int64,
                               block_size::Int64)

Compute top-k positively and negatively correlated gene pairs for ALL genes.

Strategy:
1. Full N×N dot product matrix via BLAS dgemm (multi-threaded, ~4s for 28K genes)
2. Vectorized broadcast normalization to correlations (~1s)
3. Per-gene top-k extraction via partialsortperm on column views (~15s)

Total ~20s vs ~572s in R (27× speedup). The full N×N matrix is ~6GB which fits
in memory. The block_size parameter is accepted for API compatibility but ignored.
"""
function mcview_calc_gg_mc_top_cor(daf, egc_epsilon::Float64,
                                    k::Int64, block_size::Int64)
    cache = _get_cached_egc(daf, egc_epsilon)
    n_genes = length(cache.gene_names)
    n_mc = size(cache.egc_log, 2)

    # Precompute inverse row norms
    inv_norms = Vector{Float64}(undef, n_genes)
    @inbounds for g in 1:n_genes
        inv_norms[g] = cache.row_norms[g] < 1e-15 ? 0.0 : 1.0 / cache.row_norms[g]
    end

    # Phase 1: Full N×N dot product matrix via BLAS (uses all available threads)
    cor_mat = cache.egc_log * cache.egc_log'  # n_genes × n_genes

    # Phase 2: Normalize to correlations via vectorized broadcast scaling
    # cor_mat[i,j] = dot(i,j) * inv_norms[i] * inv_norms[j]
    cor_mat .*= inv_norms           # row scale: cor_mat[i,j] *= inv_norms[i]
    cor_mat .*= inv_norms'          # col scale: cor_mat[i,j] *= inv_norms[j]
    # Zero diagonal (self-correlations)
    @inbounds for g in 1:n_genes
        cor_mat[g, g] = 0.0
    end

    # Phase 3: Extract top-k per gene using partialsortperm on column views
    max_results = n_genes * 2 * k
    result_gene1 = Vector{String}(undef, max_results)
    result_gene2 = Vector{String}(undef, max_results)
    result_cor = Vector{Float64}(undef, max_results)
    result_type = Vector{String}(undef, max_results)
    result_count = 0

    for gi in 1:n_genes
        if inv_norms[gi] == 0.0
            continue
        end

        gene1_name = cache.gene_names[gi]
        col = view(cor_mat, :, gi)  # contiguous column view (no copy)

        # Top-k positive
        pos_indices = partialsortperm(col, 1:k, rev=true)
        for idx in pos_indices
            c = @inbounds col[idx]
            if !isnan(c) && c != 0.0
                result_count += 1
                @inbounds result_gene1[result_count] = gene1_name
                @inbounds result_gene2[result_count] = cache.gene_names[idx]
                @inbounds result_cor[result_count] = c
                @inbounds result_type[result_count] = "pos"
            end
        end

        # Top-k negative
        neg_indices = partialsortperm(col, 1:k, rev=false)
        for idx in neg_indices
            c = @inbounds col[idx]
            if !isnan(c) && c != 0.0
                result_count += 1
                @inbounds result_gene1[result_count] = gene1_name
                @inbounds result_gene2[result_count] = cache.gene_names[idx]
                @inbounds result_cor[result_count] = c
                @inbounds result_type[result_count] = "neg"
            end
        end
    end

    # Trim to actual size
    return Dict(
        "gene1" => result_gene1[1:result_count],
        "gene2" => result_gene2[1:result_count],
        "cor" => result_cor[1:result_count],
        "type" => result_type[1:result_count]
    )
end


# ==============================================================================
# 5. calc_marker_genes - Top marker genes per metacell (replaces ~6.7s R path)
# ==============================================================================

"""
    mcview_calc_marker_genes(daf, genes_per_metacell::Int64,
                              minimal_max_log_fraction::Float64,
                              minimal_relative_log_fraction::Float64,
                              fold_change_reg::Float64,
                              use_abs::Bool)

Compute top-k marker genes per metacell from the EGC matrix stored in DAF.

Replicates the R `calc_marker_genes` + `select_top_fold_genes_per_metacell`:
1. log2(egc + 1e-5) transform
2. Filter genes by max log-fraction threshold
3. Compute fold-change: sweep row medians
4. Apply fold_change_reg offset and minimal_relative_log_fraction filter
5. Per-metacell top-k extraction via partialsortperm

Returns Dict with "metacell", "gene", "rank", "fp" as parallel vectors.
"""
function mcview_calc_marker_genes(daf, genes_per_metacell::Int64,
                                   minimal_max_log_fraction::Float64,
                                   minimal_relative_log_fraction::Float64,
                                   fold_change_reg::Float64,
                                   use_abs::Bool)
    # Step 1: Build log2(egc + 1e-5) matrix (n_genes × n_mc)
    egc_raw, gene_names, metacell_names = _build_egc_raw(daf)
    n_genes, n_mc = size(egc_raw)

    egc_log = Matrix{Float64}(undef, n_genes, n_mc)
    @inbounds for c in 1:n_mc
        for g in 1:n_genes
            egc_log[g, c] = log2(egc_raw[g, c] + 1e-5)
        end
    end

    # Step 2: Filter genes by max log-fraction threshold
    # Compute max per row (gene)
    max_log_fractions = Vector{Float64}(undef, n_genes)
    @inbounds for g in 1:n_genes
        mx = -Inf
        for c in 1:n_mc
            v = egc_log[g, c]
            if v > mx
                mx = v
            end
        end
        max_log_fractions[g] = mx
    end

    interesting_mask = max_log_fractions .>= minimal_max_log_fraction
    interesting_indices = findall(interesting_mask)
    n_interesting = length(interesting_indices)

    if n_interesting == 0
        return Dict(
            "metacell" => String[],
            "gene" => String[],
            "rank" => Int64[],
            "fp" => Float64[]
        )
    end

    # Step 3: Compute fold-change matrix for interesting genes
    # fold = egc_log[g, :] - median(egc_log[g, :])
    # Then apply fold_change_reg and filter by minimal_relative_log_fraction
    fold_matrix = Matrix{Float64}(undef, n_interesting, n_mc)
    sort_buf = Vector{Float64}(undef, n_mc)

    @inbounds for (row, g) in enumerate(interesting_indices)
        # Copy row for median computation
        for c in 1:n_mc
            sort_buf[c] = egc_log[g, c]
        end
        sort!(sort_buf)
        # Compute median
        if n_mc % 2 == 0
            med = (sort_buf[n_mc ÷ 2] + sort_buf[n_mc ÷ 2 + 1]) / 2.0
        else
            med = sort_buf[n_mc ÷ 2 + 1]
        end

        for c in 1:n_mc
            fp_val = egc_log[g, c] - med + fold_change_reg
            if fp_val < minimal_relative_log_fraction
                fold_matrix[row, c] = NaN  # Mark as below threshold
            else
                fold_matrix[row, c] = fp_val
            end
        end
    end

    # Step 4: Per-metacell top-k extraction
    k = min(genes_per_metacell, n_interesting)
    max_results = n_mc * k
    result_metacell = Vector{String}(undef, max_results)
    result_gene = Vector{String}(undef, max_results)
    result_rank = Vector{Int64}(undef, max_results)
    result_fp = Vector{Float64}(undef, max_results)
    result_count = 0

    sort_vals = Vector{Float64}(undef, n_interesting)

    for c in 1:n_mc
        # Build sort values for this metacell column
        @inbounds for row in 1:n_interesting
            v = fold_matrix[row, c]
            if isnan(v)
                sort_vals[row] = -Inf  # NaN values sort last
            elseif use_abs
                sort_vals[row] = abs(v)
            else
                sort_vals[row] = v
            end
        end

        # Get top-k indices
        actual_k = min(k, n_interesting)
        top_indices = partialsortperm(sort_vals, 1:actual_k, rev=true)

        mc_name = String(metacell_names[c])
        for (rank, idx) in enumerate(top_indices)
            fp_val = fold_matrix[idx, c]
            result_count += 1
            @inbounds result_metacell[result_count] = mc_name
            @inbounds result_gene[result_count] = String(gene_names[interesting_indices[idx]])
            @inbounds result_rank[result_count] = rank
            @inbounds result_fp[result_count] = isnan(fp_val) ? NaN : fp_val
        end
    end

    return Dict(
        "metacell" => result_metacell[1:result_count],
        "gene" => result_gene[1:result_count],
        "rank" => result_rank[1:result_count],
        "fp" => result_fp[1:result_count]
    )
end


# ==============================================================================
# 6. Pseudobulk computation for cell-level data
# ==============================================================================

"""
    mcview_compute_pseudobulk(cells_daf, group_field::AbstractString,
                               group_values::Vector{String},
                               cell_mask::Vector{Bool})

Compute pseudobulk (sum of UMIs per gene) for specified groups of cells.
Operates entirely in Julia on the sparse UMIs matrix — no per-gene R round-trips.

Arguments:
- cells_daf: DAF object containing cell × gene UMIs matrix
- group_field: Name of the cell-level vector to group by (e.g., "batch_set_id")
- group_values: Specific values to compute pseudobulk for (one column per value)
- cell_mask: Boolean vector (length = n_cells). TRUE = include this cell.

Returns Dict with "matrix" (n_genes × n_groups Float64), "gene_names", "group_names".
"""
function mcview_compute_pseudobulk(cells_daf, group_field::AbstractString,
                                    group_values::Vector{String},
                                    cell_mask::Vector{Bool})
    gene_names = axis_vector(cells_daf, "gene")
    n_genes = length(gene_names)
    n_groups = length(group_values)

    # Get the sparse UMIs matrix: cell × gene (CSC, columns = genes)
    umis = get_matrix(cells_daf, "cell", "gene", "UMIs")

    # Get group assignments
    group_vec = get_vector(cells_daf, "cell", group_field)

    # Build result matrix
    result = zeros(Float64, n_genes, n_groups)

    for (gi, gval) in enumerate(group_values)
        # Build combined mask: cells in this group AND passing the cell filter
        combined = falses(length(cell_mask))
        @inbounds for i in eachindex(cell_mask)
            combined[i] = cell_mask[i] && (String(group_vec[i]) == gval)
        end

        n_selected = sum(combined)
        if n_selected == 0
            continue
        end

        # Efficient sparse matrix-vector multiply: gene sums for selected cells
        # mask' * umis = (1 × cells) * (cells × genes) = 1 × genes
        mask_float = Float64.(combined)
        gene_sums = mask_float' * umis  # 1 × n_genes

        @inbounds for g in 1:n_genes
            result[g, gi] = gene_sums[g]
        end
    end

    return Dict(
        "matrix" => result,
        "gene_names" => String.(gene_names),
        "group_names" => group_values
    )
end
