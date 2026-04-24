# ==============================================================================
# Dataset Management Functions
# ==============================================================================

#' Remove a dataset
#'
#' @param dataset name of the dataset to remove
#'
#' @examples
#' \dontrun{
#' dataset_rm("PBMC163k")
#' }
#'
#' @export
dataset_rm <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data) && !is.null(mc_data[[dataset]])) {
        mc_data[[dataset]] <- NULL
        mcv_set("mc_data", mc_data)
        cli::cli_alert_success("Removed dataset '{dataset}'")
    } else {
        cli::cli_alert_warning("Dataset '{dataset}' not found")
    }
}

#' List available datasets
#'
#' @return Character vector of dataset names
#'
#' @examples
#' \dontrun{
#' dataset_ls()
#' }
#'
#' @export
dataset_ls <- function() {
    .Deprecated("dataset_names")
    dataset_names()
}

#' calculate the k top correlated and anti correlated genes for each gene
#'
#' @param egc egc matrix (normalized metacell counts per gene)
#' @param k number of top correlations per direction
#' @param egc_epsilon epsilon for log2 computation
#' @param daf_obj optional DAF object; when provided, tries Julia acceleration first
#'
#' @noRd
calc_gg_mc_top_cor <- function(egc, k = 30, egc_epsilon = 1e-5, daf_obj = NULL) {
    # Full N×N correlation matrix (BLAS-backed via tgs_cross_cor_blas)
    lfp <- log2(egc + egc_epsilon)

    cm <- tryCatch(
        .Call(
            "tgs_cross_cor_blas",
            t(lfp),
            t(lfp),
            TRUE,
            FALSE,
            FALSE,
            0,
            new.env(parent = parent.frame()),
            PACKAGE = "tgstat"
        ),
        error = function(e) {
            prev_blas <- getOption("tgs_use.blas")
            if (!isTRUE(prev_blas)) {
                options(tgs_use.blas = TRUE)
                on.exit(options(tgs_use.blas = prev_blas), add = TRUE)
            }
            tgs_cor(t(lfp), y = t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)
        }
    )

    prev_max_proc <- getOption("tgs_max.processes")
    if (is.null(prev_max_proc) || prev_max_proc != 1) {
        options(tgs_max.processes = 1)
        on.exit(options(tgs_max.processes = prev_max_proc), add = TRUE)
    }

    top_cor <- tgs_knn(cm, knn = k, diag = FALSE) %>%
        rename(gene1 = col1, gene2 = col2, cor = val) %>%
        mutate(type = "pos")
    anti_cor <- tgs_knn(-cm, knn = k, diag = FALSE) %>%
        rename(gene1 = col1, gene2 = col2, cor = val) %>%
        mutate(cor = -cor) %>%
        mutate(type = "neg")

    gg_mc_top_cor <- bind_rows(top_cor, anti_cor)

    gg_mc_top_cor <- gg_mc_top_cor %>%
        filter(!is.na(cor))

    return(gg_mc_top_cor)
}
