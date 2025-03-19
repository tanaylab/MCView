test_bundle <- function(project_bundle) {
    expect_true(fs::dir_exists(project_bundle))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "cache")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "cache", "PBMC163k")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "config")))
    expect_true(fs::file_exists(fs::path(project_bundle, "app.R")))
    expect_true(fs::file_exists(fs::path(project_bundle, "project", "config", "config.yaml")))
    expect_true(fs::file_exists(fs::path(project_bundle, "project", "config", "about.Rmd")))

    test_dataset(fs::path(project_bundle, "project"), "PBMC163k")

    config_file <- project_config_file(project_dir)
    config <- yaml::read_yaml(config_file)
    config_bundle <- yaml::read_yaml(fs::path(project_bundle, "project", "config", "config.yaml"))

    expect_equal(config, config_bundle)
}


test_that("create_bundle works", {
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC")

    project_bundle <- fs::path(bundle_dir, "PBMC")
    withr::defer(fs::dir_delete((project_bundle)))
    test_bundle(project_bundle)
})

test_that("create_bundle works with self-contained=TRUE", {
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC_sc", self_contained = TRUE)

    project_bundle <- fs::path(bundle_dir, "PBMC_sc")
    test_bundle(project_bundle)
})

test_that("create_bundle works with overwrite=TRUE", {
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC")
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC", overwrite = TRUE)

    project_bundle <- fs::path(bundle_dir, "PBMC")
    withr::defer(fs::dir_delete((project_bundle)))
    test_bundle(project_bundle)
})

test_that("create_bundle fails with overwrite=FALSE", {
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC")
    expect_error(MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC"))
})

test_that(
    "bundle app launches",
    {
        withr::local_dir(fs::path(bundle_dir, "PBMC"))

        skip_on_ci()
        skip_on_cran()
        skip_on_travis()
        skip_on_appveyor()
        x <- processx::process$new(
            "Rscript",
            "app.R"
        )
        Sys.sleep(5)
        expect_true(x$is_alive())
        x$kill()
    }
)

test_that(
    "bundle app launches when self-contaiend=TRUE",
    {
        withr::local_dir(fs::path(bundle_dir, "PBMC_sc"))

        skip_on_ci()
        skip_on_cran()
        skip_on_travis()
        skip_on_appveyor()
        x <- processx::process$new(
            "Rscript",
            "app.R"
        )
        Sys.sleep(5)
        expect_true(x$is_alive())
        x$kill()
    }
)
