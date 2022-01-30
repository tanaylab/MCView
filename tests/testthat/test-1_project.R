test_that("create_project works", {
    MCView::create_project(project_dir)

    config_file <- project_config_file(project_dir)
    expect_true(fs::file_exists(config_file))
    expect_true(fs::dir_exists(fs::path(project_dir, "config")))
    expect_true(fs::dir_exists(project_cache_dir(project_dir)))
    expect_true(fs::file_exists(project_about_file(project_dir)))
    expect_true(fs::file_exists(project_help_file(project_dir)))

    config <- yaml::read_yaml(config_file)

    required_fields <- c(
        "title",
        "tabs",
        "help"
    )

    purrr::walk(required_fields, ~ expect_true(rlang::has_name(config, .x)))
})
