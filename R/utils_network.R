

mct_propagate_flow_through_metacell <- function(mct, m) {
    time_points <- which(mct@mc_t[m, ] > 0)

    m_t_dist <- mct@mc_t[m, ]

    probs_trans_all <- list()
    probs_trans_all$probs <- matrix(0, nrow = nrow(mct@mc_t), ncol = ncol(mct@mc_t))
    probs_trans_all$step_m <- list()

    for (i in 1:(ncol(mct@mc_t) - 1)) {
        probs_trans_all$step_m[[i]] <- matrix(0, nrow = nrow(mct@mc_t), ncol = nrow(mct@mc_t))
    }

    for (t in time_points) {
        p_mc <- rep(0, nrow(mct@mc_t))
        p_mc[m] <- m_t_dist[t]

        probs_trans <- mct_propagate_backward_from_t(mct, t, p_mc)

        probs_trans_all$probs <- probs_trans_all$probs + probs_trans$probs

        for (i in 1:max(1, t - 1)) {
            probs_trans_all$step_m[[i]] <- probs_trans_all$step_m[[i]] + probs_trans$step_m[[i]]
        }
    }

    return(probs_trans_all)
}


mct_propagate_backward_from_t <- function(mct, t, mc_p) {
    max_t <- ncol(mct@mc_t)
    step_m <- list()
    for (i in 1:(max_t - 1)) {
        step_m[[i]] <- Matrix::Matrix(0, nrow = nrow(mct@mc_t), ncol = nrow(mct@mc_t), sparse = T)
    }

    probs <- matrix(0, nrow = nrow(mct@mc_t), ncol = ncol(mct@mc_t))
    probs[, t] <- mc_p
    if (t > 1) {
        for (i in (t - 1):1) {
            step_m[[i]] <- Matrix::Matrix(t(t(as.matrix(mct@mc_backward[[i]])) *
                probs[, i + 1]), sparse = T)
            probs[, i] <- as.matrix(mct@mc_backward[[i]]) %*%
                probs[, i + 1]
        }
    }

    return(list(probs = probs, step_m = step_m))
}

mctnetwork_mc_rank_from_color_ord <- function(mct_id, mc_colors, colors_ordered = NULL) {
    mct <- scdb_mctnetwork(mct_id)
    if (is.null(mct)) {
        cli_abort("cannot find mctnet object {mct_id} when trying to plot net flows")
    }
    mc <- scdb_mc(mct@mc_id)
    mc_mean_age <- apply(mct@mc_t, 1, function(x) {
        return(mean(x * (1:length(x)) / sum(x)))
    })
    col_to_rank <- c(1:length(colors_ordered))
    names(col_to_rank) <- colors_ordered

    mc_rank <- 1000 * col_to_rank[mc_colors] + mc_mean_age
    mc_rank <- rank(mc_rank)
    names(mc_rank) <- as.character(1:length(mc_rank))

    mc_flows <- mctnetwork_get_flow_mat(mct, -1)
    diag(mc_flows) <- 0
    mc_flow_fore <- mc_flows / rowSums(mc_flows)
    mc_flow_back <- t(mc_flows) / colSums(mc_flows)

    mc_rank1 <- mc_rank

    # cluster flow graph
    # rank clusters by mean rank
    # rank mc by membership + time
    if (0) {
        for (i in 1:20) {
            rank_fore <- mc_flow_fore %*% mc_rank1
            rank_back <- mc_flow_back %*% mc_rank1

            dfore <- rank_fore - mc_rank1
            dback <- rank_back - mc_rank1

            f <- sign(dfore) == sign(dback)
            f[is.na(f)] <- F
            dlt <- sign(dfore) * pmin(abs(dfore), abs(dback))
            dlt[!f] <- 0
            mc_rank1 <- rank(mc_rank1 + dlt * 0.5)
            names(mc_rank1) <- as.character(1:length(mc_rank))
        }
    }

    return(mc_rank1)
}
