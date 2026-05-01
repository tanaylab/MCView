#' MCView color palettes
#'
#' Named list of canonical palettes used by MCView plots. Centralizes
#' previously-duplicated color spec so changes propagate.
#'
#' Entries:
#' \describe{
#'   \item{`expression`}{11-step diverging RdBu (low→high), used for gene-expression
#'         heatmaps and 2D projection color-by-gene.}
#'   \item{`expression_with_zero`}{`white` followed by `viridis::viridis_pal()(6)`,
#'         used for "raw" expression with a zero-anchored low end.}
#'   \item{`cell_type_score`}{`white`→Reds→`black` ramp, used for cell-type score plots.}
#'   \item{`edge`}{Single near-black hex, used for plot edges / outlines that were
#'         previously hardcoded to `"black"`.}
#' }
#'
#' @export
mcview_palette <- list(
    expression = c(
        "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
        "#F7F7F7",
        "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
    ),
    expression_with_zero = c("white", viridis::viridis_pal()(6)),
    cell_type_score = c("white", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"),
    edge = "#1f2933"
)
