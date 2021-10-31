MCView::create_project(project_dir)

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

test_that("cell_metadata_to_metacell works", {
    set.seed(60427)
    n_cells <- 5e6
    cell_metadata <- tibble(
        cell_id = 1:n_cells,
        md1 = sample(1:5, size = n_cells, replace = TRUE),
        md2 = rnorm(n = n_cells)
    )

    cell_To_metacell <- tibble(
        cell_id = 1:n_cells,
        metacell = sample(0:1535, size = n_cells, replace = TRUE)
    )

    metadata <- cell_metadata_to_metacell(cell_metadata, cell_To_metacell)

    expect_equivalent(metadata$md1[1:5], c(
        3.01056763285024, 3.05427782888684, 2.99559748427673, 3.02988361119849,
        2.96641221374046
    ))
    expect_equivalent(metadata$md2[1:5], c(
        0.00148788798193235, 0.0181742754379464, 0.0271565244845097,
        0.0409454954704143, 0.0102880587793541
    ))
    expect_equal(metadata$metacell, 0:1535)
    expect_equal(colnames(metadata), c("metacell", "md1", "md2"))
})

test_that("cell_metadata_to_metacell works with categorical variable", {
    set.seed(60427)
    n_cells <- 5e6
    cell_metadata <- tibble(
        cell_id = 1:n_cells,
        md1 = sample(paste0("batch", 1:5), size = n_cells, replace = TRUE)
    )

    cell_To_metacell <- tibble(
        cell_id = 1:n_cells,
        metacell = sample(0:1535, size = n_cells, replace = TRUE)
    )

    metadata <- cell_metadata_to_metacell(cell_metadata, cell_To_metacell, categorical = TRUE)

    expect_equivalent(metadata$batch1[1:5], c(
        0.201421800947867, 0.213386552041756, 0.198876053699657, 0.206535141800247,
        0.201992753623188
    ))
    expect_equivalent(metadata$batch2[1:5], c(
        0.198755924170616, 0.195271722443967, 0.191695285669685, 0.200369913686806,
        0.191425120772947
    ))
    expect_equal(metadata$metacell, 0:1535)
    expect_equal(colnames(metadata), c("metacell", paste0("batch", 1:5)))
    sums <- metadata %>%
        column_to_rownames("metacell") %>%
        as.matrix() %>%
        rowSums()
    expect_equivalent(sums, rep(1, length(sums)))
})

test_that("cell_metadata_to_metacell_from_h5ad works", {
    # Skip on ci
    h5ad_file <- "/home/aviezerl/src/metacell.shiny/cells.h5ad"
    skip_if(!fs::file_exists(h5ad_file))

    metadata <- cell_metadata_to_metacell_from_h5ad(h5ad_file, c("pile", "properly_sampled_cell"))
    expect_true(all(metadata$properly_sampled_cell == 1))
    expect_equal(metadata$pile[1:5], c(0, 0, 1, 1, 1))
    expect_equal(metadata$metacell, 0:1535)
    expect_equal(colnames(metadata), c("metacell", "pile", "properly_sampled_cell"))
})

test_that("cell_metadata_to_metacell_from_h5ad works with categorical=TRUE", {
    # Skip on ci
    h5ad_file <- "/home/aviezerl/src/metacell.shiny/cells.h5ad"
    skip_if(!fs::file_exists(h5ad_file))

    metadata <- cell_metadata_to_metacell_from_h5ad(h5ad_file, c("batch"), categorical = TRUE)

    expect_setequal(
        colnames(metadata),
        c(
            "metacell", "8batchall", "b_cells", "cd14", "cd34",
            "cd4_t_helper", "cd56_nk", "cyto_t", "memory_t", "naive_cyto",
            "naive_t", "reg_t"
        )
    )

    sums <- metadata %>%
        column_to_rownames("metacell") %>%
        as.matrix() %>%
        rowSums()
    expect_equivalent(sums, rep(1, length(sums)))
    expect_equal(metadata$metacell, 0:1535)
    expect_equal(metadata$cd34[1:5], c(1, 1, 1, 1, 0))
})

# TODO:
# test edge cases cell_metadata_to_metacell


test_that("cell_metadata_to_metacell_from_h5ad works with custom function", {
    # Skip on ci
    h5ad_file <- "/home/aviezerl/src/metacell.shiny/cells.h5ad"
    skip_if(!fs::file_exists(h5ad_file))

    metadata <- cell_metadata_to_metacell_from_h5ad(h5ad_file, c("pile", "properly_sampled_cell"), func = function(x) median(x) * 2)
    expect_true(all(metadata$properly_sampled_cell == 2))
    expect_equal(metadata$pile[1:5], c(0, 0, 2, 2, 2))
    expect_equal(metadata$metacell, 0:1535)
    expect_equal(colnames(metadata), c("metacell", "pile", "properly_sampled_cell"))
})


test_that("import_dataset works with metadata fields", {
    set.seed(60427)
    dataset <- "PBMC163k_md"

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata_fields = c("grouped", "pile", "candidate")
    )
    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>%
        as_tibble() %>%
        mutate(metacell = as.character(metacell))

    library(anndata)
    adata <- anndata::read_h5ad(fs::path(raw_dir, "metacells.h5ad"))

    expect_equivalent(
        adata$obs %>%
            rownames_to_column("metacell") %>%
            select(metacell, grouped, pile, candidate),
        saved_md
    )

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("import_dataset works with metadata fields fails with non existing fields", {
    set.seed(60427)
    dataset <- "PBMC163k_md"

    expect_error(
        MCView::import_dataset(
            project = project_dir,
            dataset = dataset,
            anndata_file = fs::path(raw_dir, "metacells.h5ad"),
            cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
            metadata_fields = "savta"
        )
    )
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
        expect_warning(
            MCView::import_dataset(
                project = project_dir,
                dataset = dataset,
                anndata_file = fs::path(raw_dir, "metacells.h5ad"),
                cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
                metadata = metadata_df
            )
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

test_that("import_dataset works with metadata colors but without metadata", {
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
            metadata_colors = metadata_colors
        )
    )
})

# update metadata and colors

test_that("update_metadata works", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df %>% select(-md3)
    )

    metadata_modified <- metadata_df %>%
        mutate(md1 = rep(5, n()))

    update_metadata(
        project = project_dir,
        dataset = dataset,
        metadata = metadata_modified,
        overwrite = FALSE
    )

    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equal(saved_md$md1, metadata_modified$md1)
    expect_equivalent(saved_md, metadata_modified)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("update_metadata works with overwrite=TRUE", {
    set.seed(60427)
    dataset <- "PBMC163k_md"
    metadata_df <- get_test_metadata(raw_dir)

    MCView::import_dataset(
        project = project_dir,
        dataset = dataset,
        anndata_file = fs::path(raw_dir, "metacells.h5ad"),
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv"),
        metadata = metadata_df %>% select(-md3)
    )

    metadata_modified <- metadata_df %>%
        select(metacell) %>%
        mutate(md1 = rep(5, n()))

    update_metadata(
        project = project_dir,
        dataset = dataset,
        metadata = metadata_modified,
        overwrite = TRUE
    )

    test_dataset(project_dir, dataset)
    test_cache_file_exists("metadata.tsv", project_dir, dataset)

    saved_md <- fread(fs::path(project_dir, "cache", dataset, "metadata.tsv")) %>% as_tibble()
    expect_equivalent(saved_md, metadata_modified)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("update_metadata_colors work", {
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

    metadata_colors_modified <- metadata_colors
    metadata_colors_modified$md1 <- NULL
    metadata_colors_modified$md2 <- list(
        colors = c("red", "red", "red"),
        breaks = c(1, 2, 3)
    )

    update_metadata_colors(
        project_dir,
        dataset,
        metadata_colors_modified,
        overwrite = FALSE
    )

    saved_md_colors <- qs::qread(fs::path(project_dir, "cache", dataset, "metadata_colors.qs"))
    expect_equivalent(metadata_colors_modified, saved_md_colors)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("update_metadata_colors work with overwrite = TRUE", {
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

    metadata_colors_modified <- metadata_colors
    metadata_colors_modified$md1 <- NULL
    metadata_colors_modified$md2 <- list(
        colors = c("red", "red", "red"),
        breaks = c(1, 2, 3)
    )

    update_metadata_colors(
        project_dir,
        dataset,
        metadata_colors_modified,
        overwrite = TRUE
    )

    saved_md_colors <- qs::qread(fs::path(project_dir, "cache", dataset, "metadata_colors.qs"))
    expect_equivalent(metadata_colors_modified, saved_md_colors)

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})

test_that("update_metadata_colors fails without metadata", {
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
        cell_type_colors_file = fs::path(raw_dir, "cluster-colors.csv")
    )

    prev_metadata_file <- fs::path(project_cache_dir(project_dir), dataset, "metadata.tsv")
    if (fs::file_exists(prev_metadata_file)) {
        fs::file_delete(prev_metadata_file)
    }

    expect_error(
        update_metadata_colors(
            project_dir,
            dataset,
            metadata_colors,
            overwrite = TRUE
        )
    )

    withr::defer(unlink(fs::path(project_dir, "cache", dataset)))
    withr::defer(gc())
})
