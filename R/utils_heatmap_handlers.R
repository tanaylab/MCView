# utils_heatmap_handlers.R - Heatmap tooltip + download handlers
#
# Split from R/utils_heatmap.R (2026-05-01).  Output handlers for the
# heatmap module: hover-tooltip and download (matrix + gene list + clipboard).
# Companions: utils_heatmap_data.R, utils_heatmap_server.R,
# utils_heatmap_ui.R, utils_heatmap_help.R.

#' Build the hover tooltip UI for the heatmap
#'
#' Extracted from heatmap_reactives to reduce function complexity.
#' Renders the hover_info output that displays gene/metacell details on hover.
#'
#' @param output Shiny output object
#' @param mat Reactive expression returning the heatmap matrix
#' @param metacell_filter Reactive value for metacell filtering
#' @param metacell_types Reactive expression returning metacell type annotations
#' @param cell_type_colors Reactive expression returning cell type color mappings
#' @param dataset Reactive expression returning the current dataset
#' @param input Shiny input object
#' @param globals Reactive values for global state
#' @param mode Character string indicating heatmap mode
#' @param gene_modules Reactive expression returning gene module data
#' @param genes Reactive expression returning selected genes (can be NULL)
#' @noRd
heatmap_tooltip_handler <- function(output, mat, metacell_filter, metacell_types, cell_type_colors, dataset, input, globals, mode, gene_modules, genes, applied_params) {
    output$hover_info <- renderUI({
        m <- mat()
        req(m)
        req(input$heatmap_hover)
        req(metacell_types())
        req(cell_type_colors())

        hover <- input$heatmap_hover
        m <- filter_heatmap_by_metacell(m, metacell_filter())
        gene <- get_gene_by_heatmap_coord(m, hover$y)
        metacell <- get_metacell_by_heatmap_coord(m, hover$x)
        value <- m[gene, metacell]

        req(metacell)
        req(gene)

        # taken from https://gitlab.com/-/snippets/16220
        left_px <- hover$coords_css$x
        top_px <- hover$coords_css$y

        mcell_stats <- metacell_types() %>%
            filter(metacell == !!metacell)

        style <- glue(
            "position:absolute; z-index:100; left: {left_px + 2}px; top: {top_px + 2}px;"
        )

        top_genes <- mcell_stats %>%
            glue::glue_data("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")

        gene_prefix <- "Gene"
        if (mode == "Gene modules" && !(gene %in% genes())) {
            gene_prefix <- "Gene module"
        }

        gene_name <- gene_label(gene, dataset(), gene_modules())
        if (mode == "Gene modules") {
            gene_name <- gene
        }

        mcell_tooltip <- paste(
            glue("{gene_prefix}: {gene_name}"),
            glue("Metacell: {metacell}"),
            glue("Value: {round(value, digits=2)}"),
            glue("Cell type: {mcell_stats$cell_type}"),
            glue("Top genes: {top_genes}"),
            ifelse(has_name(mcell_stats, "mc_age"), glue("Metacell age (E[t]): {round(mcell_stats$mc_age, digits=2)}"), ""),
            sep = "<br/>"
        )

        metadata <- get_markers_metadata(dataset, applied_params()$selected_md, metacell_types, globals)
        if (!is.null(metadata)) {
            mc_md <- metadata %>%
                filter(metacell == !!metacell) %>%
                select(-metacell)
            md_tooltip <- purrr::map_chr(colnames(mc_md), ~
                glue("{.x}: {mc_md[[.x]][1]}")) %>%
                paste(collapse = "<br/>")
            mcell_tooltip <- paste(mcell_tooltip, md_tooltip)
        }

        wellPanel(
            style = style,
            p(HTML(mcell_tooltip))
        )
    })
}


#' Set up download handlers for the heatmap
#'
#' Extracted from heatmap_reactives to reduce function complexity.
#' Creates download handlers for the marker matrix and gene list.
#'
#' @param output Shiny output object
#' @param mat Reactive expression returning the heatmap matrix
#' @param markers Reactive value containing the current marker genes
#' @param metacell_filter Reactive value for metacell filtering
#' @param dataset Reactive expression returning the current dataset
#' @param input Shiny input object
#' @param metacell_types Reactive expression returning metacell type annotations
#' @param globals Reactive values for global state
#' @noRd
heatmap_download_handlers <- function(output, mat, markers, metacell_filter, dataset, input, metacell_types, globals, applied_params, ns) {
    output$download_matrix <- downloadHandler(
        filename = function() {
            paste("markers_matrix-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            m <- mat()[rev(1:nrow(mat())), , drop = FALSE]
            if (length(metacell_filter()) > 0) {
                m <- m[, intersect(colnames(m), metacell_filter()), drop = FALSE]
            }
            if (input$include_metadata) {
                metadata <- get_markers_metadata(dataset, applied_params()$selected_md, metacell_types, globals)
                metadata_m <- metadata %>%
                    as.data.frame() %>%
                    column_to_rownames("metacell") %>%
                    t() %>%
                    as.matrix()
                ct_m <- metacell_types() %>%
                    select(metacell, cell_type) %>%
                    deframe()
                ct_m <- ct_m[colnames(m)]
                metadata_m <- rbind(t(as.matrix(ct_m)), metadata_m)
                rownames(metadata_m)[1] <- "Cell type"
                m <- rbind(m, metadata_m)
            }
            fwrite(
                m %>%
                    as.data.frame() %>%
                    rownames_to_column("gene"),
                file,
                row.names = FALSE
            )
        }
    )

    output$download_genes <- downloadHandler(
        filename = function() {
            paste("markers-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                data.frame(gene = markers()),
                file,
                row.names = FALSE,
                col.names = FALSE
            )
        }
    )

    # Render copy genes button using utility function
    output$copy_genes_button <- clipboard_copy_button_ui(
        ns, "copy_genes_to_clipboard", markers,
        label = "Copy genes to clipboard",
        style = "margin: 5px 5px 5px 15px; background-color: #17a2b8; color: white; border: none;",
        tooltip = "Copy all genes to system clipboard"
    )

    clipboard_copy_button_server(
        input, "copy_genes_to_clipboard", markers, globals,
        message_template = "Copied {count} genes to clipboard"
    )
}
