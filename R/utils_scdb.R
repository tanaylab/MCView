init_temp_scdb <- function(config, dataset, scdb_dir = fs::path(tempdir(), dataset, "scdb"), force_copy = TRUE) {
    scdb_dir <<- scdb_dir
    fs::dir_create(scdb_dir)
    scdb_init(base_dir = scdb_dir, force_reinit = TRUE)

    if (!force_copy && dir.exists(system.file("scdb", package = "MCView"))) {
        scdb_files <- list.files(system.file("scdb", package = "MCView"), full.names = TRUE)
    } else {
        all_files <- list.files(config$scdb)
        metacells <- unique(c(config$matrix, config$mc, config$mc2d, config$network))
        files_to_copy <- purrr::map(metacells, ~ grep(glue(".+\\.{.x}\\.Rda"), all_files, value = TRUE)) %>% do.call(c, .)
        scdb_files <- fs::path(config$scdb, files_to_copy)
    }

    purrr::walk(scdb_files, ~ fs::link_create(.x, fs::path(scdb_dir, basename(.x))))
}
