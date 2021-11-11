
test_that("update_metacell_types works", {
    metacell_types <- fread(fs::path(raw_dir, "metacell-types.csv"))
    metacell_types$cell_type[1:50] <- 1
    test_file <- fs::path(raw_dir, "metacell-types-test.csv")
    fwrite(metacell_types, test_file)
    withr::defer(unlink(test_file))
    MCView::update_metacell_types(project_dir, "PBMC163k", test_file)

    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    metacell_types <- get_metacell_types_data("PBMC163k")

    expect_true(all(metacell_types$cell_type[1:50] == 1))
})

test_that("update_metacell_types works with missing cell types", {
    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    orig_mc_types <- get_metacell_types_data("PBMC163k")
    metacell_types <- fread(fs::path(raw_dir, "metacell-types.csv"))
    metacell_types <- metacell_types %>% filter(!(metacell %in% 1:50))
    test_file <- fs::path(raw_dir, "metacell-types-test.csv")
    fwrite(metacell_types, test_file)
    withr::defer(unlink(test_file))
    MCView::update_metacell_types(project_dir, "PBMC163k", test_file)

    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    metacell_types <- get_metacell_types_data("PBMC163k")
    f <- metacell_types$metacell %in% 1:50
    expect_equal(sum(f), 50)
    expect_true(all(metacell_types$cell_type[f] == "(Missing)"))
    expect_true(all(!is.na(metacell_types$cell_type[!f])))

    expect_equal(nrow(metacell_types %>% filter(!f) %>% anti_join(metacell_types %>% filter(!f))), 0)
})

# TODO: Test for cell types that do not exist in colors

test_that("update_cell_type_colors works", {
    cell_type_colors <- fread(fs::path(raw_dir, "cluster-colors.csv"))
    cell_type_colors$color[1:5] <- "black"
    test_file <- fs::path(raw_dir, "cluster-colors-test.csv")
    fwrite(cell_type_colors, test_file)
    withr::defer(unlink(test_file))
    MCView::update_cell_type_colors(project_dir, "PBMC163k", test_file)

    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    cell_type_colors <- get_cell_type_data("PBMC163k")

    expect_true(all(cell_type_colors$color[1:5] == "black"))
})

# TODO: Test for cell types that do not exist in metacell types
