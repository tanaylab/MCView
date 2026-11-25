# The names of the genes and metacells of mc_mat are stored in mc_dimnames, so that the app gets them without
# reading the matrix.

local_mc_data <- function(env = parent.frame()) {
    prev <- mget("mc_data", envir = globalenv(), ifnotfound = list(NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())), envir = env)
}

fake_mc_mat <- function(genes = c("g1", "g2", "g3"), metacells = c("M1", "M2")) {
    matrix(seq_len(length(genes) * length(metacells)), length(genes), length(metacells), dimnames = list(genes, metacells))
}

# a dataset directory with mc_mat and mc_sum, without mc_dimnames (as imported before it existed)
legacy_dataset <- function(cache_dir, dataset = "ds", mc_mat = fake_mc_mat(), sub = NULL) {
    dir <- fs::path(cache_dir, dataset, sub %||% "")
    fs::dir_create(dir)
    qs2::qs_save(mc_mat, fs::path(dir, "mc_mat.qs2"))
    qs2::qs_save(as.numeric(colSums(mc_mat)), fs::path(dir, "mc_sum.qs2"))
    dir
}

# make reading mc_mat fail, keeping its time (a newer mc_mat makes mc_dimnames stale)
break_mc_mat <- function(dir) {
    file <- fs::path(dir, "mc_mat.qs2")
    mtime <- file.mtime(file)
    writeLines("not a qs2 file", file)
    Sys.setFileTime(file, mtime)
}

test_that("serialize_shiny_data writes the names of mc_mat next to it", {
    cache_dir <- withr::local_tempdir()
    serialize_shiny_data(fake_mc_mat(), "mc_mat", dataset = "ds", cache_dir = cache_dir)
    expect_identical(qs2::qs_read(fs::path(cache_dir, "ds", "mc_dimnames.qs2")), list(genes = c("g1", "g2", "g3"), metacells = c("M1", "M2")))
})

test_that("gene and metacell names and mc_sum use mc_dimnames without reading mc_mat", {
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir)
    serialize_shiny_data(fake_mc_mat(), "mc_mat", dataset = "ds", cache_dir = cache_dir)
    break_mc_mat(dir)

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "g3"))
    expect_identical(dataset_metacell_names("ds"), c("M1", "M2"))
    expect_identical(names(get_mc_sum("ds")), c("M1", "M2"))
    expect_error(get_mc_data("ds", "mc_mat")) # it was never read
    # a dataset or atlas without mc_mat
    expect_null(gene_names("ds", atlas = TRUE))
    expect_null(gene_names("no_such_dataset"))
})

test_that("without mc_dimnames, the names come from mc_mat and are written for the next processes", {
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir)
    atlas_dir <- legacy_dataset(cache_dir, mc_mat = fake_mc_mat(genes = c("a1", "a2"), metacells = "A1"), sub = "atlas")

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "g3"))
    expect_identical(names(get_mc_sum("ds")), c("M1", "M2"))
    expect_identical(gene_names("ds", atlas = TRUE), c("a1", "a2"))
    expect_identical(dataset_metacell_names("ds", atlas = TRUE), "A1")
    expect_identical(qs2::qs_read(fs::path(dir, "mc_dimnames.qs2")), list(genes = c("g1", "g2", "g3"), metacells = c("M1", "M2")))
    expect_identical(qs2::qs_read(fs::path(atlas_dir, "mc_dimnames.qs2")), list(genes = c("a1", "a2"), metacells = "A1"))
    # no temporary file left
    expect_identical(sort(list.files(dir, all.files = TRUE, no.. = TRUE)), c("atlas", "mc_dimnames.qs2", "mc_mat.qs2", "mc_sum.qs2"))

    # the next process doesn't read mc_mat
    break_mc_mat(dir)
    break_mc_mat(atlas_dir)
    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "g3"))
    expect_identical(gene_names("ds", atlas = TRUE), c("a1", "a2"))
})

test_that("mc_dimnames is not written to a read-only directory, silently, and the names are kept in memory", {
    skip_on_os("windows")
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir)
    Sys.chmod(dir, "555", use_umask = FALSE)
    withr::defer(Sys.chmod(dir, "755", use_umask = FALSE))
    skip_if(file.access(dir, 2) == 0, "can write to a read-only directory (root)")

    load_all_data(cache_dir, datasets = "ds")
    expect_silent(genes <- gene_names("ds"))
    expect_identical(genes, c("g1", "g2", "g3"))
    expect_false(fs::file_exists(fs::path(dir, "mc_dimnames.qs2")))
    expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), c("mc_mat.qs2", "mc_sum.qs2"))
    expect_identical(get_mc_data("ds", "mc_dimnames"), list(genes = c("g1", "g2", "g3"), metacells = c("M1", "M2")))
    expect_false(write_mc_dimnames(list(genes = "x", metacells = "y"), dir))
    expect_false(write_mc_dimnames(list(genes = "x", metacells = "y"), fs::path(dir, "no_such_dir")))
})

test_that("mc_dimnames is not written to a dataset linked from outside the project", {
    skip_on_os("windows")
    local_mc_data()
    elsewhere <- withr::local_tempdir()
    dir <- legacy_dataset(elsewhere, dataset = "atlas_ds")
    cache_dir <- withr::local_tempdir()
    fs::link_create(dir, fs::path(cache_dir, "ds"))

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "g3"))
    expect_false(fs::file_exists(fs::path(dir, "mc_dimnames.qs2")))
})

test_that("stale mc_dimnames is corrected once mc_mat is read", {
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir, mc_mat = fake_mc_mat(genes = c("g1", "g2", "new")))
    qs2::qs_save(list(genes = c("g1", "g2", "old"), metacells = c("M1", "M2")), fs::path(dir, "mc_dimnames.qs2"))

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "old")) # mc_mat wasn't read
    mc_mat <- get_mc_data("ds", "mc_mat")
    expect_identical(gene_names("ds"), c("g1", "g2", "new"))
    expect_identical(qs2::qs_read(fs::path(dir, "mc_dimnames.qs2"))$genes, c("g1", "g2", "new"))

    # the same names are not written again
    mtime <- file.mtime(fs::path(dir, "mc_dimnames.qs2"))
    Sys.setFileTime(fs::path(dir, "mc_dimnames.qs2"), mtime + 100)
    load_all_data(cache_dir, datasets = "ds")
    mc_mat <- get_mc_data("ds", "mc_mat")
    expect_equal(as.numeric(file.mtime(fs::path(dir, "mc_dimnames.qs2"))), as.numeric(mtime + 100))
})

test_that("mc_dimnames older than mc_mat is ignored", {
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir, mc_mat = fake_mc_mat(genes = c("g1", "g2", "new")))
    qs2::qs_save(list(genes = c("g1", "g2", "old"), metacells = c("M1", "M2")), fs::path(dir, "mc_dimnames.qs2"))
    # a new mc_mat copied over a cache that has mc_dimnames
    Sys.setFileTime(fs::path(dir, "mc_dimnames.qs2"), file.mtime(fs::path(dir, "mc_mat.qs2")) - 100)

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "new"))
    expect_identical(qs2::qs_read(fs::path(dir, "mc_dimnames.qs2"))$genes, c("g1", "g2", "new"))
})

test_that("an unreadable mc_dimnames falls back to mc_mat", {
    local_mc_data()
    cache_dir <- withr::local_tempdir()
    dir <- legacy_dataset(cache_dir)
    writeLines("not a qs2 file", fs::path(dir, "mc_dimnames.qs2"))

    load_all_data(cache_dir, datasets = "ds")
    expect_identical(gene_names("ds"), c("g1", "g2", "g3"))
    expect_identical(qs2::qs_read(fs::path(dir, "mc_dimnames.qs2"))$metacells, c("M1", "M2"))
})

test_that("writing mc_dimnames doesn't change the key of the shiny cache", {
    prev <- mget(c("config", "project"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer({
        purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv()))
        shiny::shinyOptions(cache = NULL)
    })
    root <- withr::local_tempdir()
    withr::local_envvar(MCVIEW_SHINY_CACHE_DIR = root)
    p <- withr::local_tempdir()
    fs::dir_create(fs::path(p, "config"))
    writeLines("title: x", fs::path(p, "config", "config.yaml"))
    legacy_dataset(fs::path(p, "cache"))
    assign("project", p, envir = globalenv())
    assign("config", list(), envir = globalenv())

    config_shiny_cache()
    expect_length(list.dirs(root, recursive = FALSE), 1)
    write_mc_dimnames(list(genes = "x", metacells = "y"), fs::path(p, "cache", "ds"))
    config_shiny_cache()
    expect_length(list.dirs(root, recursive = FALSE), 1)
})
