

get_tab_steps <- function(tab, module_id = NULL) {
    steps <- tibble(element = names(help_config[[tab]]), intro = help_config[[tab]])

    if (!is.null(module_id)) {
        steps <- steps %>% mutate(element = NS(module_id)(element))
    }

    steps <- steps %>%
        mutate(element = ifelse(element == "NA", NA, paste0("#", element)))

    steps <- as.data.frame(steps)

    return(steps)
}
