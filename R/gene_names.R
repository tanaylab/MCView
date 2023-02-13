modify_gene_names <- function(orig_names, gene_names = NULL) {
    if (is.null(gene_names)) {
        return(orig_names)
    }
    purrr::walk(c("gene_name", "alt_name"), ~ {
        if (!has_name(gene_names, .x)) {
            cli_abort("{.field gene_names} doesn't have the required field {.field {.x}}")
        }
    })

    res <- tibble(gene_name = orig_names) %>%
        left_join(gene_names, by = "gene_name") %>%
        mutate(new_name = ifelse(!is.na(alt_name), alt_name, gene_name)) %>%
        pull(new_name)

    return(res)
}
