get_test_metadata <- function(raw_dir) {
    metacell_types <- fread(fs::path(raw_dir, "metacell-types.csv")) %>% as_tibble()

    metadata <- metacell_types %>%
        select(-cell_type) %>%
        mutate(
            md1 = rnorm(1:n()),
            md2 = sample(1:n(), n()),
            md3 = md2 / sum(md2),
            md4 = rep(1, n())
        )

    return(metadata)
}
