# Adapted from https://stefanengineering.com/2019/07/06/delete-rows-from-shiny-dt-datatable/

#' A column of delete buttons for each row in the data frame for the first column
#'
#' @param df data frame
#' @return A DT::datatable with escaping turned off that has the delete buttons in the first column and \code{df} in the other
#'
#' @noRd
delete_button_column <- function(df, id, ns, ...) {
    # function to create one action button as string
    f <- function(i) {
        # https://shiny.rstudio.com/articles/communicating-with-js.html
        as.character(actionButton(paste(id, i, sep = "_"),
            label = NULL, icon = icon("trash"),
            onclick = paste0('Shiny.setInputValue(\"', ns("delete_pressed"), '\",  this.id, {priority: "event"})')
        ))
    }

    delete_col <- unlist(lapply(df$metacell, f))

    # Return a data table
    DT::datatable(cbind(df, delete = delete_col),
        # Need to disable escaping for html as string to work
        escape = FALSE,
        options = list(
            columnDefs = list(list(targets = 2, sortable = FALSE)),
            pageLength = 20
        )
    )
}

dt_selector_column <- function(df, id, ns, choices, ...) {
    # function to create selector
    f <- function(i, selected) {
        selectInput(ns(paste(id, i, sep = "_")), choices = choices, selected = selected, width = "150px", label = "") %>% as.character()
    }

    select_col <- purrr::map2(df$metacell, df$cell_type, f) %>% unlist()

    return(cbind(df %>% select(-cell_type), select = select_col))
}


#' Extracts the row id number from the id string
#' @param idstr the id string formated as id_INDEX
#'
#' @noRd
#' @return INDEX from the id string id_INDEX
parse_delete_event <- function(idstr) {
    res <- as.integer(sub(".*_([0-9]+)", "\\1", idstr))
    if (!is.na(res)) res
}
