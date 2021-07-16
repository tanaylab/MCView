raw_dir <- test_path("raw")

fs::dir_create(raw_dir)

download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", fs::path(raw_dir, "PBMC_processed.tar.gz"))
untar(fs::path(raw_dir, "PBMC_processed.tar.gz"), exdir = raw_dir)

project_dir <- test_path("PBMC")
bundle_dir <- test_path("bundle")

withr::defer(
    {
        unlink(raw_dir)
        unlink(project_dir)
        unlink(bundle_dir)
    },
    teardown_env()
)

required_files <- c("mc_mat.qs", "mc_sum.qs", "mc2d.qs", "gg_mc_top_cor.qs", "metacell_types.tsv", "cell_type_colors.tsv")

test_dataset <- function(project_dir, dataset) {
    dataset_dir <- fs::path(project_dir, "cache", dataset)
    expect_true(fs::dir_exists(project_dir))
    expect_true(fs::dir_exists(dataset_dir))
    purrr::walk(required_files, ~ expect_true(fs::file_exists(fs::path(dataset_dir, .x))))
}
