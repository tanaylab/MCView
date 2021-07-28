MCView::create_project("PBMC", project_dir = project_dir)

test_that("import_dataset works with metadata dataframe", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata file", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)
    metadata_file <- fs::path(test_dir, "metadata.csv")
    tgutil::fwrite(metadata_df, metadata_file)

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_file
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata without all metacells", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_df <- metadata_df[-(1:10), ]

    expect_warning(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df
        )
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset fails with metadata with unknown metacells", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir) %>% mutate(metacell = as.character(metacell))
    metadata_df <- bind_rows(
        metadata_df %>% mutate(metacell = as.character(metacell)),
        tibble(metacell = "savta", md1 = 1, md2 = 2, md3 = 3, md4 = 2)
    )
    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df
        )
    )
})

test_that("import_dataset fails when metadata is categorical", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir) %>%
        mutate(md1 = "savta")

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df
        )
    )
})


test_that("import_dataset works with metadata colors", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md1 = list(
            colors = c("darblue", "white", "darkred"),
            breaks = c(0, 0.5, 1)
        ),
        md2 = list(
            colors = c("black", "gray", "yellow"),
            breaks = c(5, 500, 1000)
        ),
        md3 = list(
            colors = c("green", "purple", "orange"),
            breaks = c(0, 0.5, 1)
        )
    )

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df,
        metadata_colors = metadata_colors
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    saved_md_colors <- qs::qread(fs::path(project_dir, "cache", dataset, "metadata_colors.qs"))
    expect_equivalent(metadata_colors, saved_md_colors)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata colors file", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md1 = list(
            colors = c("darblue", "white", "darkred"),
            breaks = c(0, 0.5, 1)
        ),
        md2 = list(
            colors = c("black", "gray", "yellow"),
            breaks = c(5, 500, 1000)
        ),
        md3 = list(
            colors = c("green", "purple", "orange"),
            breaks = c(0, 0.5, 1)
        )
    )

    metadata_colors_file <- fs::path(test_dir, "metadata_colors.yaml")
    yaml::write_yaml(metadata_colors, metadata_colors_file)

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df,
        metadata_colors = metadata_colors_file
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    saved_md_colors <- qs::qread(fs::path(project_dir, "cache", dataset, "metadata_colors.qs"))
    expect_equivalent(metadata_colors, saved_md_colors)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata colors without breaks", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md1 = list(
            colors = c("darblue", "white", "darkred")
        ),
        md2 = list(
            colors = c("black", "gray", "yellow"),
            breaks = c(5, 500, 1000)
        ),
        md3 = list(
            colors = c("green", "purple", "orange")
        )
    )

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df,
        metadata_colors = metadata_colors
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(metadata_df, saved_md)

    saved_md_colors <- qs::qread(fs::path(project_dir, "cache", dataset, "metadata_colors.qs"))
    expect_equivalent(metadata_colors, saved_md_colors)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata colors with wrong number of breaks", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md1 = list(
            colors = c("darblue", "white", "darkred")
        ),
        md2 = list(
            colors = c("black", "gray", "yellow"),
            breaks = c(5, 500, 1000, 5000, 8000)
        ),
        md3 = list(
            colors = c("green", "purple", "orange")
        )
    )

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df,
            metadata_colors = metadata_colors
        )
    )
})

test_that("import_dataset works with metadata colors with wrong number of colors", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md1 = list(
            colors = c("darblue", "white", "darkred")
        ),
        md2 = list(
            colors = c("black", "gray", "yellow", "red"),
            breaks = c(5, 500, 1000)
        ),
        md3 = list(
            colors = c("green", "purple", "orange")
        )
    )

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df,
            metadata_colors = metadata_colors
        )
    )
})


test_that("import_dataset works with metadata colors with non-existing metadata fields", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    metadata_colors <- list(
        md18 = list(
            colors = c("darblue", "white", "darkred")
        )
    )

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata = metadata_df,
            metadata_colors = metadata_colors
        )
    )
})
