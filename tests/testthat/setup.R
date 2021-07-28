set.seed(60427)

test_dir <- tempdir()

raw_dir <- fs::path(test_dir, "raw")

fs::dir_create(raw_dir)

# set the timeout to more than 60 seconds in order to be able to download the files on github
options(timeout = 1e4)
download.file(
    "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz",
    fs::path(raw_dir, "PBMC_processed.tar.gz")
)

untar(fs::path(raw_dir, "PBMC_processed.tar.gz"), exdir = raw_dir)

project_dir <- fs::path(test_dir, "PBMC")
bundle_dir <- fs::path(test_dir, "bundle")

withr::defer(
    {
        unlink(raw_dir)
        unlink(project_dir)
        unlink(bundle_dir)
        unlink(test_dir)
    },
    teardown_env()
)

required_files <- c(
    "mc_mat.qs",
    "mc_sum.qs",
    "mc2d.qs",
    "gg_mc_top_cor.qs",
    "metacell_types.tsv",
    "cell_type_colors.tsv"
)

test_cache_file_exists <- function(file, project_dir, dataset) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::file_exists(fs::path(dataset_dir, file)))
}

test_dataset <- function(project_dir, dataset) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::dir_exists(project_dir))
    expect_true(fs::dir_exists(dataset_dir))
    purrr::walk(required_files, ~ test_cache_file_exists(.x, project_dir, dataset))
}
