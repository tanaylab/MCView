test_that("import_dataset works", {
    set.seed(60427)
    dataset <- "PBMC163k"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = NULL
    )
    test_dataset(project_dir, dataset)
})

test_that("dataset_ls works", {
    expect_equal(dataset_ls(project_dir), c("PBMC163k"))
})

test_that("import_dataset works with implicit cell type field", {
    set.seed(60427)
    dataset <- "PBMC163k_ctf"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad")
    )
    test_dataset(project_dir, dataset)
    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with non existing cell type field", {
    set.seed(60427)
    dataset <- "PBMC163k_ctf1"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = "savta"
    )
    test_dataset(project_dir, dataset)
    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with cell type field and cell_type_colors file", {
    set.seed(60427)
    dataset <- "PBMC163k_ctf_ctcf"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )
    test_dataset(project_dir, dataset)
    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metacell_types file", {
    set.seed(60427)
    dataset <- "PBMC163k_mct"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        metacell_types_file = fs::path(raw_dir, "metacell-types.csv")
    )
    test_dataset(project_dir, dataset)
    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metacell_types file and cell_type_colors file", {
    set.seed(60427)
    dataset <- "PBMC163k_mct"
    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        metacell_types_file = fs::path(raw_dir, "metacell-types.csv"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )
    test_dataset(project_dir, dataset)
    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("dataset_rm works", {
    dataset_rm(project_dir, "PBMC163k_ctf")
    expect_true(!fs::dir_exists(fs::path(project_dir, "cache", "PBMC163k_ctf")))
    expect_true(!("PBMC163k_ctf" %in% dataset_ls(project_dir)))
})
