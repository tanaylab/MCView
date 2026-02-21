dt_selector_column <- function(df, id, ns, choices, ...) {
    # function to create selector
    f <- function(i, selected) {
        selectInput(ns(paste(id, i, sep = "_")), choices = choices, selected = selected, width = "150px", label = "") %>% as.character()
    }

    select_col <- purrr::map2(df$metacell, df$cell_type, f) %>% unlist()

    return(cbind(df %>% select(-cell_type), select = select_col))
}
