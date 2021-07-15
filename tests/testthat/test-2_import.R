required_files <- c("mc_mat.qs", "mc_sum.qs", "mc2d.qs", "gg_mc_top_cor.qs", "metacell_types.tsv", "cell_type_colors.tsv")

MCView::create_project(project_dir)

test_dataset <- function(project_dir, dataset) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::dir_exists(project_dir))
    expect_true(fs::dir_exists(dataset_dir))
    purrr::walk(required_files, ~ expect_true(fs::file_exists(fs::path(dataset_dir, .x))))
}

test_that("import_dataset works", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k",
        anndata_file = fs::path(raw_dir, "metacells.h5ad")
    )
    test_dataset(project_dir, "PBMC163k")
})

test_that("import_dataset works with cell type field", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_ctf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = "cluster"
    )
    test_dataset(project_dir, "PBMC163k_ctf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf")))
})

test_that("import_dataset works with cell type field and cell_type_colors file", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_ctf_ctcf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = "cluster",
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )
    test_dataset(project_dir, "PBMC163k_ctf_ctcf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf_ctcf")))
})

test_that("import_dataset works with metacell_types file", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_mct",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        metacell_types_file = fs::path(raw_dir, "metacell-types.csv")
    )
    test_dataset(project_dir, "PBMC163k_mct")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_mct")))
})

test_that("import_dataset works with metacell_types file and cell_type_colors file", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_mct_ctcf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        metacell_types_file = fs::path(raw_dir, "metacell-types.csv"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )
    test_dataset(project_dir, "PBMC163k_mct_ctcf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_mct_ctcf")))
})
