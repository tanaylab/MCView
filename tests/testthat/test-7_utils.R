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
