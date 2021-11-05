init_temp_scdb <- function(scdb, matrix, mc, mc2d, network = NULL, dataset, scdb_dir = fs::path(tempdir(), dataset, "scdb"), force_copy = TRUE) {
    scdb_dir <<- scdb_dir
    fs::dir_create(scdb_dir)
    metacell::scdb_init(base_dir = scdb_dir, force_reinit = TRUE)

    scdb <- normalizePath(scdb)

    if (!force_copy && dir.exists(system.file("scdb", package = "MCView"))) {
        scdb_files <- list.files(system.file("scdb", package = "MCView"), full.names = TRUE)
    } else {
        all_files <- list.files(scdb)
        metacells <- unique(c(matrix, mc, mc2d, network))
        files_to_copy <- purrr::map(metacells, ~ grep(glue(".+\\.{.x}\\.Rda"), all_files, value = TRUE)) %>% do.call(c, .)
        scdb_files <- fs::path(scdb, files_to_copy)
    }

    purrr::walk(scdb_files, ~ {
        link <- fs::path(scdb_dir, basename(.x))
        if (!fs::file_exists(link)) {
            fs::link_create(.x, link)
        }
    })
}
