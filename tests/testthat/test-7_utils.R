test_that("subset_dist matches subsetting the full distance matrix", {
    x <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("m", 1:20), NULL))
    d <- dist(x)
    expect_identical(subset_dist(d, labels(d)), d)
    for (labs in list(paste0("m", c(2, 5, 6, 11, 20)), paste0("m", c(7, 3, 9)), paste0("m", 1:19))) {
        expect_equal(as.matrix(subset_dist(d, labs)), as.matrix(d)[labs, labs])
    }
})

test_that("cor_dist matches tgs_dist(tgs_cor())", {
    set.seed(60427)
    feat <- matrix(rnorm(30 * 700), nrow = 30, dimnames = list(NULL, paste0("m", 1:700)))
    expected <- tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE))
    for (blas in c(TRUE, FALSE)) { # with tgstat's BLAS cor_dist computes the distances differently
        withr::local_options(tgs_use.blas = blas)
        d <- cor_dist(feat)
        expect_equal(as.vector(d), as.vector(expected))
        expect_identical(attributes(d)[c("Size", "Labels", "Diag", "Upper", "method", "class")], attributes(expected)[c("Size", "Labels", "Diag", "Upper", "method", "class")])
        expect_identical(fastcluster::hclust(d, method = "ward.D2")$merge, fastcluster::hclust(expected, method = "ward.D2")$merge)
    }

    feat[, 5] <- 1 # constant metacell: falls back to tgs_dist(tgs_cor())
    expect_identical(cor_dist(feat), tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)))
})

test_that("cluster_metacells returns the same clustering with and without the app cache", {
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    hc <- cluster_metacells(feat, "key")
    shiny::shinyOptions(cache = cachem::cache_mem())
    withr::defer(shiny::shinyOptions(cache = NULL))
    expect_identical(cluster_metacells(feat, "key")$order, hc$order)
    expect_identical(cluster_metacells(feat, "key")$order, hc$order) # cache hit
})

test_that("cluster_metacells uses a stored clustering without reading the distances", {
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    stored <- list(key = fastcluster::hclust(dist(t(feat))))
    expect_identical(cluster_metacells(feat, "key", cached_dist = function() stop("read the distances"), precomputed_clusterings = stored), stored$key)
    d <- dist(t(feat))
    expect_identical(cluster_metacells(feat, "other", cached_dist = function() d, precomputed_clusterings = stored)$merge, fastcluster::hclust(d, method = "ward.D2")$merge)
})

test_that("cluster_metacells shares clusterings through MCVIEW_CLUSTERING_CACHE_DIR", {
    withr::local_envvar(MCVIEW_CLUSTERING_CACHE_DIR = withr::local_tempdir())
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    hc <- cluster_metacells(feat, "12345")
    expect_length(dir(Sys.getenv("MCVIEW_CLUSTERING_CACHE_DIR")), 1)
    # the same key in another process (no app cache, nothing stored with the project): read, not computed
    expect_identical(cluster_metacells(feat, "12345", cached_dist = function() stop("computed it again")), hc)
})

test_that("calc_diff_expr p-values match chisq.test", {
    set.seed(60427)
    # counts and projected (non-integer) counts, with many genes far from the diagonal
    mat <- cbind(a = rpois(3000, rexp(3000, 0.01)), b = rexp(3000, 0.01) * runif(3000, 0, 5))
    rownames(mat) <- paste0("g", 1:3000)
    egc <- t(t(mat) / colSums(mat)) + 1e-5
    df <- calc_diff_expr(mat, egc, c("a", "b"))
    f <- abs(df$diff) >= 1.5
    expect_gt(sum(f), 500)
    m <- mat[df$gene[f], ]
    tots <- colSums(m)
    expected <- apply(m, 1, function(x) suppressWarnings(stats::chisq.test(matrix(c(x, tots), nrow = 2))$p.value))
    expect_identical(df$pval[f], unname(expected))
    expect_true(all(is.na(df$pval[!f])))

    # a single gene far from the diagonal
    mat1 <- rbind(g1 = c(a = 500, b = 3), g2 = c(a = 5000, b = 300))
    df1 <- calc_diff_expr(mat1, t(t(mat1) / colSums(mat1)) + 1e-5, c("a", "b"))
    expect_identical(df1$pval[df1$gene == "g1"], suppressWarnings(stats::chisq.test(matrix(c(500, 3, 500, 3), nrow = 2))$p.value))
})

test_that("calc_obs_exp_mc_df tests the counts of the chosen metacell", {
    set.seed(60427)
    obs <- matrix(rpois(2000 * 6, rexp(2000 * 6, 0.01)), 2000, 6, dimnames = list(paste0("g", 1:2000), paste0("mc", 1:6)))
    proj <- obs * runif(length(obs), 0, 6)
    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    assign("mc_data", list(toy = list2env(list(mc_mat = obs, mc_sum = colSums(obs), projected_mat = proj, projected_mat_sum = colSums(proj)))), envir = globalenv()) # an environment per dataset, as in the app
    assign("egc_epsilon", 1e-5, envir = globalenv())

    df <- calc_obs_exp_mc_df("toy", "mc4")
    f <- !is.na(df$pval)
    expect_gt(sum(f), 100)
    m <- cbind(obs[df$gene[f], "mc4"], proj[df$gene[f], "mc4"])
    tots <- colSums(m)
    expected <- apply(m, 1, function(x) suppressWarnings(stats::chisq.test(matrix(c(x, tots), nrow = 2))$p.value))
    expect_equal(df$pval[f], unname(expected))
    expect_equal(df$Observed, unname(obs[df$gene, "mc4"] / sum(obs[, "mc4"]) + 1e-5))
})

test_that("the diff. expression hover shows the diff rounded to 3 digits", {
    prev <- mget("expr_breaks", envir = globalenv(), ifnotfound = list(NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    assign("expr_breaks", c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1), envir = globalenv())
    df <- tibble(gene = c("a", "b", "c"), gene_name = gene, m1 = c(1.25e-5, 0.01, 0.002), m2 = c(0.0149, 0.001, 0.002), diff = c(-10.22512, 3.32193, 0), pval = c(0, 1e-3, NA), col = factor(c("darkblue", "darkred", "gray"), levels = c("gray", "darkred", "darkblue")))
    p <- plot_mc_mc_gene(df, "m1", "m2")
    # set_hovertemplate writes the as.character() of each field
    expect_identical(as.character(p$data$diff_text[match(df$gene, p$data$gene)]), c("-10.225", "3.322", "0"))
})

test_that("set_hovertemplate keeps each point's id first in its customdata", {
    df <- data.frame(id = c("a", "b", "c"), v = c(1.5, NA, 3))
    fig <- plotly::plot_ly(df, x = 1:3, y = 1:3, customdata = ~ rev(id), type = "scatter", mode = "markers")
    tr <- set_hovertemplate(fig, df, "id", "v", "%{customdata[0]}: %{customdata[1]}")$x$data[[1]]
    expect_identical(unname(tr$customdata), rbind(c("c", "3"), c("b", "NA"), c("a", "1.5")))
    expect_identical(tr$hovertemplate, "%{customdata[0]}: %{customdata[1]}")
    expect_null(tr$text)

    # fields with one value for all the trace's points go in its template
    df$p <- "Not computed"
    tr <- set_hovertemplate(fig, df, "id", c("p", "v"), "%{customdata[0]}: %{customdata[2]} (%{customdata[1]})")$x$data[[1]]
    expect_identical(unname(tr$customdata), rbind(c("c", "3"), c("b", "NA"), c("a", "1.5")))
    expect_identical(tr$hovertemplate, "%{customdata[0]}: %{customdata[1]} (Not computed)")

    # unless plotly would read them as a placeholder
    df$p <- "50%{y}"
    tr <- set_hovertemplate(fig, df, "id", c("p", "v"), "%{customdata[0]}: %{customdata[2]} (%{customdata[1]})")$x$data[[1]]
    expect_identical(unname(tr$customdata), rbind(c("c", "50%{y}", "3"), c("b", "50%{y}", "NA"), c("a", "50%{y}", "1.5")))
    expect_identical(tr$hovertemplate, "%{customdata[0]}: %{customdata[2]} (%{customdata[1]})")

    # a point that isn't in df keeps its id
    fig <- plotly::plot_ly(x = 1:2, y = 1:2, customdata = c("a", "z"), type = "scatter", mode = "markers")
    tr <- set_hovertemplate(fig, df, "id", "v", "%{customdata[0]}: %{customdata[1]}")$x$data[[1]]
    expect_identical(unname(tr$customdata), rbind(c("a", "1.5"), c("z", "NA")))
})

test_that("plotly_event_ids reads the ids of clicked and selected points", {
    event <- function(json) jsonlite::parse_json(json, simplifyVector = TRUE) # as plotly::event_data()
    expect_identical(plotly_event_ids(event('[{"pointNumber":3,"customdata":["M1","a","b"]}]')), "M1")
    expect_identical(plotly_event_ids(event('[{"pointNumber":3,"customdata":"M1"},{"pointNumber":4,"customdata":"M2"}]')), c("M1", "M2"))
    expect_identical(plotly_event_ids(event('[{"pointNumber":3,"customdata":["M1","a"]},{"pointNumber":0,"customdata":"M2"},{"pointNumber":1}]')), c("M1", "M2", NA))
})

test_that("lazy_load_files reads the .qs2 when a .qs of the same name exists", {
    cache_dir <- withr::local_tempdir()
    fs::dir_create(fs::path(cache_dir, "ds"))
    qs2::qs_save(1:3, fs::path(cache_dir, "ds", "x.qs2"))
    writeLines("not a qs file", fs::path(cache_dir, "ds", "x.qs"))
    env <- new.env()
    lazy_load_files(env, c("x.qs", "x.qs2"), "ds", cache_dir)
    expect_identical(ls(env), "x")
    expect_identical(env$x, 1:3)
})

test_that("densify_if_mostly_nonzero densifies a dgCMatrix above 2/3 non-zeros", {
    m <- Matrix::Matrix(c(1, 2, 3, 5, 0, 4), 2, 3, sparse = TRUE)
    expect_identical(densify_if_mostly_nonzero(m), as.matrix(m))
    m <- Matrix::Matrix(c(1, 2, 3, 0, 0, 4), 2, 3, sparse = TRUE) # exactly 2/3
    expect_identical(densify_if_mostly_nonzero(m), m)
})

test_that("update_default_clusterings skips only the datasets linked from outside the project", {
    project <- withr::local_tempdir()
    fs::dir_create(fs::path(project, c("config", "cache/inside")))
    fs::file_create(fs::path(project, c("config/config.yaml", "cache/inside/default_clusterings.qs2")))
    out <- paste(cli::cli_fmt(update_default_clusterings(project)), collapse = " ")
    expect_match(out, "inside already has default_clusterings")
    expect_no_match(out, "symlink")

    skip_on_os("windows") # creating a symlink needs privileges there
    fs::link_create(withr::local_tempdir(), fs::path(project, "cache", "linked"))
    out <- paste(cli::cli_fmt(update_default_clusterings(project)), collapse = " ")
    expect_match(out, "inside already has default_clusterings")
    expect_match(out, "linked is a symlink to .* outside the project")
})

test_that("caches imported with MCView 0.2.7 to 0.2.9 are read in the current format", {
    prev <- mget("mc_data", envir = globalenv(), ifnotfound = list(NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    cache_dir <- withr::local_tempdir()
    ds_dir <- fs::path(cache_dir, "old")
    fs::dir_create(ds_dir)
    mcs <- c("M1", "M2", "M3")
    # metacells as rownames, and no max inner-fold without lateral genes
    qs2::qs_save(data.frame(umis = c(10, 20, 30), cells = 1:3, max_inner_fold = c(5, 3, 1), row.names = mcs), fs::path(ds_dir, "mc_qc_metadata.qs2"))
    qs2::qs_save(Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 1, 2, 3), x = c(5, 2, 3, 1), dims = c(3, 3), dimnames = list(c("lat", "g1", "g2"), mcs)), fs::path(ds_dir, "inner_fold_mat.qs2"))
    qs2::qs_save("lat", fs::path(ds_dir, "lateral_genes.qs2"))
    # `lateral` instead of `type`
    qs2::qs_save(tibble(gene = c("lat", "g1"), zero_fold = c(2, 1), avg = -5, obs = 3, exp = 1, metacell = "M1", lateral = c(TRUE, FALSE)), fs::path(ds_dir, "gene_zero_fold.qs2"))
    data.table::fwrite(data.frame(metacell = mcs, age = 1:3), fs::path(ds_dir, "metadata.tsv"), sep = "\t")
    # the first zero-fold format: no gene expression and no gene type
    older_dir <- fs::path(cache_dir, "older")
    fs::dir_create(older_dir)
    qs2::qs_save(matrix(c(2, 2, 0, 1, 1, 2, 0, 3, 1), 3, 3, dimnames = list(c("lat", "g1", "g2"), mcs)), fs::path(older_dir, "mc_mat.qs2"))
    qs2::qs_save(c(4, 4, 4), fs::path(older_dir, "mc_sum.qs2"))
    qs2::qs_save("lat", fs::path(older_dir, "lateral_genes.qs2"))
    qs2::qs_save(tibble(gene = c("lat", "g1"), zero_fold = c(2, 1), obs = 3, exp = 1, metacell = "M1"), fs::path(older_dir, "gene_zero_fold.qs2"))

    load_all_data(cache_dir, datasets = c("old", "older"))
    qc <- get_mc_data("old", "mc_qc_metadata")
    expect_identical(qc$metacell, mcs)
    expect_equal(qc$max_inner_fold_no_lateral, c(2, 3, 0))
    expect_equal(get_metadata("old")$umis, c(10, 20, 30))
    expect_identical(get_mc_data("old", "gene_zero_fold")$type, c("lateral", "other"))
    zero_fold <- get_mc_data("older", "gene_zero_fold")
    expect_identical(zero_fold$type, c("lateral", "other"))
    expect_equal(zero_fold$avg, log2(c(0.25, 0.5) + 1e-5))
})

test_that("MCVIEW_SHINY_CACHE_DIR keeps a shiny cache per project version across processes", {
    prev <- mget(c("config", "project"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer({
        purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv()))
        shiny::shinyOptions(cache = NULL)
    })
    root <- withr::local_tempdir()
    withr::local_envvar(MCVIEW_SHINY_CACHE_DIR = root)
    projects <- fs::path(withr::local_tempdir(), c("a", "b"))
    for (p in projects) {
        fs::dir_create(fs::path(p, c("config", "cache/ds")))
        writeLines("title: x", fs::path(p, "config", "config.yaml"))
        writeLines("metacell\tx", fs::path(p, "cache", "ds", "mc2d.tsv"))
    }
    start <- function(p) { # what run_app() does for the cache
        assign("project", p, envir = globalenv())
        assign("config", list(), envir = globalenv())
        config_shiny_cache()
        shiny::getShinyOption("cache")
    }
    n_dirs <- function() length(list.dirs(root, recursive = FALSE))

    cache <- start(projects[1])
    cache$set("plot", list(value = 1, error = FALSE))
    cache$set("failed", list(value = NULL, error = TRUE))
    expect_equal(n_dirs(), 1)
    # the next process of the same project
    cache <- start(projects[1])
    expect_equal(n_dirs(), 1)
    expect_equal(cache$get("plot")$value, 1)
    expect_true(cachem::is.key_missing(cache$get("failed"))) # errors are kept in memory only
    # another project with the same files
    start(projects[2])
    expect_equal(n_dirs(), 2)
    # a changed file
    writeLines("metacell\tx\ty", fs::path(projects[1], "cache", "ds", "mc2d.tsv"))
    cache <- start(projects[1])
    expect_equal(n_dirs(), 3)
    expect_true(cachem::is.key_missing(cache$get("plot")))
    # the config's shiny_cache_dir comes first
    assign("config", list(shiny_cache_dir = FALSE), envir = globalenv())
    config_shiny_cache()
    expect_equal(n_dirs(), 3)
})

test_that("the shiny cache directories not used for 30 days are removed", {
    prev <- mget(c("config", "project"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer({
        purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv()))
        shiny::shinyOptions(cache = NULL)
    })
    root <- withr::local_tempdir()
    old <- fs::path(root, c("shiny_cache_old_1", "other"))
    recent <- fs::path(root, "shiny_cache_recent_1")
    fs::dir_create(c(old, recent))
    writeLines("x", fs::path(old[1], "x.qs2"))
    Sys.setFileTime(old, Sys.time() - 31 * 24 * 3600)
    Sys.setFileTime(recent, Sys.time() - 29 * 24 * 3600)
    p <- withr::local_tempdir()
    fs::dir_create(fs::path(p, c("config", "cache")))
    assign("project", p, envir = globalenv())
    assign("config", list(shiny_cache_dir = root), envir = globalenv())
    config_shiny_cache()
    expect_false(fs::dir_exists(old[1]))
    expect_true(fs::dir_exists(old[2])) # not a shiny cache directory
    expect_true(fs::dir_exists(recent))
    expect_length(list.dirs(root, recursive = FALSE), 3) # and this project's
})

test_that("available_tabs leaves out the tabs a dataset has no data for", {
    prev <- mget(c("mc_data", "tab_defs"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    init_tab_defs()
    ds <- new.env()
    ds$inner_fold_mat <- matrix(1)
    assign("mc_data", list(ds = ds), envir = globalenv())
    tabs <- available_tabs("ds")
    expect_true(all(c("QC", "Manifold", "Inner-fold") %in% tabs))
    expect_false(any(c("Atlas", "Query", "Projected-fold", "Samples", "Stdev-fold", "Outliers", "Flow") %in% tabs))
    ds$atlas <- new.env()
    expect_true(all(c("Atlas", "Query", "Projected-fold") %in% available_tabs("ds")))
})

test_that("order_mc_by_most_var_genes orders all the metacells when some have too few non-zero genes", {
    set.seed(60427)
    feat <- matrix(abs(rnorm(8 * 30)) + 0.1, nrow = 8, dimnames = list(paste0("g", 1:8), paste0("m", 1:30)))
    feat[, c(4, 17)] <- 0 # metacells without non-zero genes are clustered separately
    feat[1, 4] <- 1
    ord <- order_mc_by_most_var_genes(feat, log_transform = FALSE)
    expect_setequal(ord, 1:30)
    expect_setequal(tail(ord, 2), c(4, 17))
})

test_that("initial point sizes take a few values", {
    sizes <- c(
        purrr::map_dbl(c(1280, 1440, 1920, 2560, 3840), ~ initial_proj_point_size_value(588, .x, .x * 9 / 16, weight = 0.6)),
        purrr::map_dbl(c(1280, 1440, 1920, 2560, 3840), ~ initial_proj_point_size_value(5320, .x, .x * 9 / 16))
    )
    expect_true(all(sizes * 4 == round(sizes * 4)))
    expect_true(all(sizes >= 1 & sizes <= 3))
})
