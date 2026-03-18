set.seed(60427)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)

test_dir <- tempdir()

raw_dir <- fs::path(test_dir, "raw")

fs::dir_create(raw_dir)

# Cache test data locally to avoid re-downloading on every test run
options(timeout = 1e4)
local_cache <- fs::path(Sys.getenv("HOME"), ".cache", "mcview_test_data")
cached_tarball <- fs::path(local_cache, "PBMC_processed.tar.gz")
if (!fs::file_exists(cached_tarball)) {
    fs::dir_create(local_cache)
    download.file(
        "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz",
        cached_tarball
    )
}
fs::file_copy(cached_tarball, fs::path(raw_dir, "PBMC_processed.tar.gz"))

untar(fs::path(raw_dir, "PBMC_processed.tar.gz"), exdir = raw_dir)

project_dir <- fs::path(test_dir, "PBMC")
bundle_dir <- fs::path(test_dir, "bundle")

withr::defer(
    {
        unlink(raw_dir)
        unlink(project_dir)
        unlink(bundle_dir)
        unlink(test_dir)
        unlink(fs::path(test_dir, "PBMC1"))
    },
    teardown_env()
)

test_cache_file_exists <- function(file, project_dir, dataset) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::file_exists(fs::path(dataset_dir, file)))
}

test_dataset <- function(project_dir,
                         dataset,
                         required_files = c(
                             "mc_mat.qs",
                             "mc_sum.qs",
                             "mc2d.qs",
                             "metacell_types.tsv",
                             "cell_type_colors.tsv",
                             "marker_genes.qs"
                         )) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::dir_exists(project_dir))
    expect_true(fs::dir_exists(dataset_dir))
    purrr::walk(required_files, ~ test_cache_file_exists(.x, project_dir, dataset))
}
