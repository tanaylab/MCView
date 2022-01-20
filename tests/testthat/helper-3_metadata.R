get_test_metadata <- function(raw_dir) {
    metacell_types <- tgutil::fread(fs::path(raw_dir, "metacell-types.csv")) %>% tibble::as_tibble()

    metadata <- metacell_types %>%
        dplyr::select(-cell_type) %>%
        dplyr::mutate(
            md1 = rnorm(1:n()),
            md2 = sample(1:n(), n()),
            md3 = md2 / sum(md2),
            md4 = rep(1, n()),
            md5 = sample(c("type1", "type2", "type3"), n(), replace = TRUE)
        )

    return(metadata)
}
