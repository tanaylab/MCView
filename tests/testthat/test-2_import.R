MCView::create_project("PBMC", project_dir = project_dir)

test_that("import_dataset works", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = NULL
    )
    test_dataset(project_dir, "PBMC163k")
})

test_that("dataset_ls works", {
    expect_equal(dataset_ls(project_dir), c("PBMC163k"))
})

test_that("import_dataset works with implicit cell type field", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_ctf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad")
    )
    test_dataset(project_dir, "PBMC163k_ctf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf")))
    withr::defer(gc())
})

test_that("import_dataset works with non existing cell type field", {
    set.seed(60427)
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_ctf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_field = "savta"
    )
    test_dataset(project_dir, "PBMC163k_ctf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf")))
    withr::defer(gc())
})

test_that("import_dataset works with cell type field and cell_type_colors file", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_ctf_ctcf",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )
    test_dataset(project_dir, "PBMC163k_ctf_ctcf")
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf_ctcf")))
    withr::defer(gc())
})

test_that("import_dataset works with calc_gg_cor=FALSE", {
    MCView::import_dataset(
        project = project_dir,
        dataset = "PBMC163k_gg_cor",
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        calc_gg_cor = FALSE
    )
    test_dataset(
        project_dir,
        "PBMC163k_gg_cor",
        required_files = c(
            "mc_mat.qs",
            "mc_sum.qs",
            "mc2d.qs",
            "metacell_types.tsv",
            "cell_type_colors.tsv"
        )
    )
    withr::defer(unlink(fs::path(project_dir, "cache", "PBMC163k_ctf_ctcf")))
    withr::defer(gc())
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
    withr::defer(gc())
})

test_that("import_dataset fails with non existing metacells in metacell types", {
    mc_types_file <- tempfile()
    mc_types <- tgutil::fread(fs::path(raw_dir, "metacell-types.csv"))
    mc_types$metacell[1] <- "savta"
    tgutil::fwrite(mc_types, mc_types_file)

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = "PBMC163k_mct",
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            metacell_types_file = mc_types_file
        )
    )

    withr::defer(unlink(mc_types_file))
    withr::defer(gc())
})

test_that("import_dataset warns about missing metacells in metacell types", {
    mc_types_file <- tempfile()
    mc_types <- tgutil::fread(fs::path(raw_dir, "metacell-types.csv"))
    mc_types <- mc_types[2:nrow(mc_types), ]
    tgutil::fwrite(mc_types, mc_types_file)

    expect_warning(
        MCView::import_dataset(
            project = project_dir,
            dataset = "PBMC163k_mct",
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            metacell_types_file = mc_types_file
        )
    )

    withr::defer(unlink(mc_types_file))
    withr::defer(gc())
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
    withr::defer(gc())
})

test_that("dataset_rm works", {
    dataset_rm(project_dir, "PBMC163k_ctf")
    expect_true(!fs::dir_exists(fs::path(project_dir, "cache", "PBMC163k_ctf")))
    expect_true(!("PBMC_163k_ctf" %in% dataset_ls(project_dir)))
})
