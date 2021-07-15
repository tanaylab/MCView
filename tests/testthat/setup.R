raw_dir <- test_path("raw")

fs::dir_create(raw_dir)

download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", fs::path(raw_dir, "PBMC_processed.tar.gz"))
untar(fs::path(raw_dir, "PBMC_processed.tar.gz"), exdir = raw_dir)

project_dir <- test_path("PBMC")

withr::defer(
    {
        unlink(raw_dir)
        unlink(project_dir)
    },
    teardown_env()
)
