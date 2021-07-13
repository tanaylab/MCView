mctnetwork_g_t_types <- function(net, min_time, max_time, g, mc_egc, mc_colors) {
    gene_e <- log2(1e-5 + mc_egc[g, ])

    all_types <- unique(mc_colors)

    src_mats <- list()
    targ_mats <- list()

    for (t in min_time:(max_time - 1)) {
        net_t <- net %>%
            filter(time1 == t, time2 == t + 1, type1 != "growth", type2 != "growth")


        net_t$mc_t1 <- mc_colors[as.numeric(net_t$mc1)]
        net_t$mc_t2 <- mc_colors[as.numeric(net_t$mc2)]
        net_t$src_g <- gene_e[net_t$mc1]
        net_t$targ_g <- gene_e[net_t$mc2]


        src_flow <- as.data.frame(summarize(group_by(net_t, mc_t1, mc_t2),
            g_src = weighted.mean(src_g, flow)
        ))
        targ_flow <- as.data.frame(summarize(group_by(net_t, mc_t1, mc_t2),
            g_targ = weighted.mean(targ_g, flow)
        ))


        src_mat <- as.data.frame(pivot_wider(
            data = src_flow,
            names_from = mc_t2,
            values_from = g_src,
            values_fill = list(g_src = NA)
        ))

        targ_mat <- as.data.frame(pivot_wider(
            data = targ_flow,
            names_from = mc_t2,
            values_from = g_targ,
            values_fill = list(g_targ = NA)
        ))

        rownames(src_mat) <- src_mat$mc_t1
        rownames(targ_mat) <- targ_mat$mc_t1
        src_mat <- src_mat[, -1]
        targ_mat <- targ_mat[, -1]

        src_mat <- src_mat[all_types, all_types]
        src_mat <- as.matrix(src_mat)
        src_mats[[t]] <- src_mat

        targ_mat <- targ_mat[all_types, all_types]
        targ_mat <- as.matrix(targ_mat)
        targ_mats[[t]] <- targ_mat
    }

    return(list(src_g_tt = src_mats, targ_g_tt = targ_mats))
}


add_alpha <- function(col, alpha) {
    return(rgb(t(col2rgb(col)) / 256, alpha = alpha))
}


get_sig_edge <- function(x1, x2, x2t, y1, y2, y2t, flow, col1, col2, col_alpha = 0.8) {
    x1 <- x1
    y1 <- y1
    dx <- x2t - x1
    dy <- y2t - y1

    y1t <- y1 + flow
    dxt <- x2 - x1
    dyt <- y2 - y1t

    col1 <- col2rgb(col1)[, 1]
    names(col1) <- c("red", "green", "blue")
    col2 <- col2rgb(col2)[, 1]
    names(col2) <- c("red", "green", "blue")
    res <- 0.05

    beta0 <- plogis(0, loc = 0.5, scale = 0.2)
    beta_f <- plogis(1, loc = 0.5, scale = 0.2) - plogis(0, loc = 0.5, scale = 0.2)
    polygons <- list()
    for (r in seq(0, 0.98, res)) {
        beta <- (plogis(r, loc = 0.5, scale = 0.2) - beta0) / beta_f
        beta5 <- (plogis(r + res, loc = 0.5, scale = 0.2) - beta0) / beta_f

        sx1 <- x1 + r * dx
        sy1 <- y1 + beta * dy
        sx2 <- x1 + (r + res) * dx
        sy2 <- y1 + beta5 * dy

        sx1t <- x1 + r * dxt
        sy1t <- y1t + beta * dyt
        sx2t <- x1 + (r + res) * dxt
        sy2t <- y1t + beta5 * dyt

        r_col <- r
        rgb_r <- col2["red"] * r_col + col1["red"] * (1 - r_col)
        rgb_g <- col2["green"] * r_col + col1["green"] * (1 - r_col)
        rgb_b <- col2["blue"] * r_col + col1["blue"] * (1 - r_col)
        col <- rgb(rgb_r / 256, rgb_g / 256, rgb_b / 256, col_alpha)
        poly_list <- list(
            mat = matrix(c(sx1, sx2, sx2t, sx1t, sy1, sy2, sy2t, sy1t), nrow = 2, byrow = TRUE),
            col = col,
            border = NA,
            lwd = 1
        )
        polygons <- append(polygons, list(poly_list))
    }

    return(polygons)
}


plot_vein <- function(dataset,
                      gene = NULL,
                      foc_type = NULL,
                      max_eg = -11,
                      min_eg = -15,
                      gene_shades = colorRampPalette(c("gray95", "lightblue", "blue", "darkblue", "gold"))(100),
                      min_flow = 0,
                      vein_lwd = 1,
                      color_order = NULL,
                      T_minflow_for_type = 0.005,
                      mc_annot = get_mc_data(dataset, mc_annot),
                      mc_color_key = get_mc_color_key(dataset)) {
    type_ag <- get_mc_data(dataset, "type_ag")
    type_flow <- get_mc_data(dataset, "type_flow")
    time_annot <- get_mc_data(dataset, "time_annot")
    mc_egc <- get_mc_egc(dataset)
    mc_colors <- mc_annot$mc_col
    mc_network <- get_mc_data(dataset, "mc_network")


    min_t <- min(as.numeric(colnames(type_ag)))
    max_t <- max(as.numeric(colnames(type_ag)))
    t1 <- min_t
    t2 <- max_t - 1

    if (!is.null(foc_type)) {
        key <- mc_color_key
        rownames(key) <- key$group
        foc_color <- key[foc_type, "color"]
        adj_sums <- lapply(type_flow, function(x) x[foc_color, ]) %>%
            do.call("rbind", .) %>%
            colSums()
        adj_cols <- which(adj_sums > T_minflow_for_type) %>% names()

        all_foc_colors <- c(foc_color)

        # if foc_type is terminal - plot progenitors
        if (length(adj_cols) == 1) {
            adj_sums <- lapply(type_flow, function(x) x[, foc_color]) %>%
                do.call("rbind", .) %>%
                colSums()
            adj_cols <- which(adj_sums > T_minflow_for_type) %>% names()

            all_foc_colors <- c(adj_cols)
        }
    } else {
        adj_cols <- colnames(type_flow[[t1]])
        all_foc_colors <- c(adj_cols)
    }

    if (!is.null(color_order)) {
        color_ord <- 1:length(color_order)
        names(color_ord) <- toupper(color_order)

        adj_cols <- intersect(adj_cols, color_order)
        all_foc_colors <- intersect(all_foc_colors, color_order)
        adj_cols1 <- adj_cols[order(color_ord[toupper(adj_cols)])]
        adj_cols <- adj_cols1
    }

    type_agn <- t(t(type_ag) / colSums(type_ag))
    rownames(type_agn) <- get_cell_type_colors(dataset)[rownames(type_agn)]
    foc_agn <- type_agn[adj_cols, t1:t2]

    if (!is.null(gene)) {
        a <- mctnetwork_g_t_types(mc_network, min_t, max_t, gene, mc_egc, mc_colors)
        src_g_tt <- a$src_g_tt
        targ_g_tt <- a$targ_g_tt
        plot_gene <- TRUE
        egc_to_col <- function(x) {
            if (is.null(x) | !is.finite(x)) {
                return("white")
            } else {
                xbin <- max(floor((x - min_eg) * 100 / (max_eg - min_eg)) + 1, 1)
                return(gene_shades[min(xbin, 100)])
            }
        }
    } else {
        plot_gene <- FALSE
    }

    base_y <- c(1)

    # computing the y coordinate of each vein
    for (i in 2:length(adj_cols)) {
        base_y[i] <- base_y[i - 1] + 0.2 + max(foc_agn[adj_cols[i], ] + foc_agn[adj_cols[i - 1], ])
    }
    names(base_y) <- adj_cols

    smoo_y <- list()

    get_polygons <- function() {
        polygons <- list()

        # plotting the veins, using loess smoothing on the frequencies
        for (c in adj_cols) {
            x <- t1:t2
            y <- foc_agn[c, t1:t2]
            ys <- approx(x, y, seq(t1, t2, l = 1 + (t2 - t1) * 100))
            lo <- loess(ys$y ~ ys$x, span = 0.3)$fitted
            base <- base_y[c]
            smoo_y[[c]] <- lo
            names(smoo_y[[c]]) <- ys$x
            if (plot_gene) {
                col_alpha <- 0.8
                i <- 1
                for (t in t1:(t2 - 1)) {
                    if (t == t1) {
                        col1 <- egc_to_col(src_g_tt[[t]][c, c])
                    } else {
                        col1 <- egc_to_col((targ_g_tt[[t - 1]][c, c] + src_g_tt[[t]][c, c]) / 2)
                    }
                    col2 <- egc_to_col((targ_g_tt[[t]][c, c] + src_g_tt[[t + 1]][c, c]) / 2)
                    col1 <- col2rgb(col1)[, 1]
                    names(col1) <- c("red", "green", "blue")
                    col2 <- col2rgb(col2)[, 1]
                    names(col2) <- c("red", "green", "blue")
                    for (lamda in seq(0, 0.99, l = 100)) {
                        rgb_r <- col2["red"] * lamda + col1["red"] * (1 - lamda)
                        rgb_g <- col2["green"] * lamda + col1["green"] * (1 - lamda)
                        rgb_b <- col2["blue"] * lamda + col1["blue"] * (1 - lamda)
                        col <- rgb(rgb_r / 256, rgb_g / 256, rgb_b / 256, col_alpha)

                        poly_list <- list(
                            mat = matrix(c(ys$x[i], ys$x[i + 1], ys$x[i + 1], ys$x[i], base + lo[i], base + lo[i + 1], base - lo[i + 1], base - lo[i]), nrow = 2, byrow = TRUE),
                            col = col,
                            border = NA,
                            lwd = 1
                        )
                        polygons <- append(polygons, list(poly_list))
                        i <- i + 1
                    }
                }

                poly_list <- list(
                    mat = matrix(c(ys$x, rev(ys$x), base + c(lo, rev(-lo))), nrow = 2, byrow = TRUE),
                    col = NA,
                    border = c,
                    lwd = vein_lwd
                )
                polygons <- append(polygons, list(poly_list))
            } else {
                if (!is.null(foc_type)) {
                    calpha <- add_alpha(c, 0.8)
                    poly_list <- list(
                        mat = matrix(c(ys$x, rev(ys$x), base + c(lo, rev(-lo))), nrow = 2, byrow = TRUE),
                        col = ifelse(c == foc_color, c, calpha),
                        border = NA,
                        lwd = 1
                    )
                    polygons <- append(polygons, list(poly_list))
                } else {
                    poly_list <- list(
                        mat = matrix(c(ys$x, rev(ys$x), base + c(lo, rev(-lo))), nrow = 2, byrow = TRUE),
                        col = c,
                        border = NA,
                        lwd = 1
                    )
                    polygons <- append(polygons, list(poly_list))
                }
            }
        }

        for (foc_color in all_foc_colors) {
            foc_i <- which(adj_cols == foc_color)
            base_foc <- base_y[foc_color]
            for (t in t1:(t2 - 1)) {
                flow <- type_flow[[t]]
                max_i <- length(adj_cols)
                cum_y <- -smoo_y[[foc_color]][as.character(t)]
                if (foc_i != 1) {
                    for (i in 1:(foc_i - 1)) {
                        col_i <- adj_cols[i]
                        fl <- flow[foc_color, col_i] * 2
                        if (!is.null(fl) & length(fl) > 0 & fl > min_flow) {
                            if (plot_gene) {
                                col1 <- egc_to_col(src_g_tt[[t]][foc_color, col_i])
                                col2 <- egc_to_col(targ_g_tt[[t]][foc_color, col_i])
                            } else {
                                col1 <- foc_color
                                col2 <- col_i
                            }
                            poly_list <- get_sig_edge(
                                x1 = t, x2 = t + 1, x2t = t + 1 - 2 * fl - 0.05,
                                y1 = base_foc + cum_y,
                                y2 = base_y[col_i] + smoo_y[[col_i]][as.character(t + 1)],
                                y2t = base_y[col_i] + smoo_y[[col_i]][as.character(round(t + 1 - 2 * fl - 0.05, 2))],
                                flow = fl,
                                col1 = col1, col2 = col2
                            )
                            polygons <- append(polygons, poly_list)
                            cum_y <- cum_y + fl
                        }
                    }
                }
                cum_y <- smoo_y[[foc_color]][as.character(t)]
                if (foc_i != max_i) {
                    for (i in max_i:(foc_i + 1)) {
                        col_i <- adj_cols[i]
                        fl <- flow[foc_color, col_i] * 2
                        if (fl > 0) {
                            if (plot_gene) {
                                col1 <- egc_to_col(src_g_tt[[t]][foc_color, col_i])
                                col2 <- egc_to_col(targ_g_tt[[t]][foc_color, col_i])
                            } else {
                                col1 <- foc_color
                                col2 <- col_i
                            }
                            poly_list <- get_sig_edge(
                                x1 = t, x2t = t + 1, x2 = t + 1 - 2 * fl - 0.05,
                                y1 = base_foc + cum_y - fl,
                                y2t = base_y[col_i] - smoo_y[[col_i]][as.character(t + 1)],
                                y2 = base_y[col_i] - smoo_y[[col_i]][as.character(round(t + 1 - 2 * fl - 0.05, 2))],
                                flow = fl,
                                col1 = col1, col2 = col2
                            )
                            polygons <- append(polygons, poly_list)
                            cum_y <- cum_y - fl
                        }
                    }
                }
            }
        }

        return(polygons)
    }

    # Plot

    if (plot_gene) {
        par(mar = rep(1, 4))
        layout(mat = matrix(c(1, 2), nrow = 1), widths = c(10, 50))
        plot_color_bar(c(1, 100), c(min_eg, max_eg), gene_shades, title = "Gene expression")
    }


    plot(0,
        xlim = c(t1 - 0.5, t2),
        ylim = c(min(base_y) - 0.5, max(base_y) + 0.5),
        yaxt = "n",
        xaxt = "n",
        bg = "transparent",
        xlab = "Time (E[t])",
        ylab = "Cell fraction",
        main = gene
    )

    axis(side = 1, at = t1:t2, labels = time_annot$time_desc[t1:t2])

    x_pos <- 0.5
    y_c <- mean(base_y)

    segments(x0 = x_pos, x1 = x_pos, y0 = y_c - 0.25, y1 = y_c + 0.25, lwd = 4)
    segments(x0 = x_pos - 0.05, x1 = x_pos + 0.05, y0 = y_c - 0.25, y1 = y_c - 0.25, lwd = 4)
    segments(x0 = x_pos - 0.05, x1 = x_pos + 0.05, y0 = y_c + 0.25, y1 = y_c + 0.25, lwd = 4)

    text(x = x_pos - 0.2, y = y_c, labels = "25%", cex = 1)

    plot_vein_polygons <- function(polygons) {
        # Plot polygons
        xs <- purrr::map(polygons, ~ c(.x$mat[1, ], NA)) %>%
            purrr::flatten() %>%
            purrr::map_dbl(~.x)
        ys <- purrr::map(polygons, ~ c(.x$mat[2, ], NA)) %>%
            purrr::flatten() %>%
            purrr::map_dbl(~.x)
        cols <- purrr::map_chr(polygons, ~ .x$col)
        lwds <- purrr::map_dbl(polygons, ~ .x$lwd)
        borders <- purrr::map_chr(polygons, ~ .x$border)
        polygon(xs, ys, col = cols, border = borders, lwd = lwds)
    }

    future_promise({
        get_polygons()
    }) %...>% plot_vein_polygons()
}
