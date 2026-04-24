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
    lfp <- log2(egc + egc_epsilon)
    streaming_top_k_cor(
        X = lfp, k = k, diag = FALSE, min_cor = 0
    )
}
