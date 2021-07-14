plot_propagation_net_metacell <- function(dataset,
                                          metacell,
                                          edge_w_scale = 5e-4,
                                          fr_scale = 2,
                                          max_lwd = 10,
                                          metacell_type= get_mc_data(dataset, "mc_annot")) {
    metacell <- as.character(metacell)

    mct_probs_trans <- get_mc_data(dataset, "mct_probs_trans")
    mc_network <- get_mc_data(dataset, "mc_network")
    mc_rank <- get_mc_data(dataset, "mc_rank")
    time_annot <- get_mc_data(dataset, "time_annot")

    mc_propagate <- mct_probs_trans[[metacell]]$step_m

    dx_back <- 0
    dy_ext <- 0
    bg_col <- "gray90"
    bg <- "white"

    score_shades <- colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000)

    nn <- mc_network %>%
        filter(flow > 1e-4) %>%
        mutate(
            x1 = time1,
            x2 = time2,
            y1 = as.numeric(mc_rank[as.character(mc1)]),
            y2 = as.numeric(mc_rank[as.character(mc2)]),
            x1 = ifelse(type1 == "growth", x1 + 0.3, x1),
            x2 = ifelse(type2 == "growth", x2 + 0.3, x2),
            x1 = ifelse(type1 == "norm_b" | type1 == "extend_b", x1 - dx_back, x1),
            x2 = ifelse(type2 == "norm_b" | type2 == "extend_b", x2 - dx_back, x2),
            y1 = ifelse(type1 == "src", max(y1) / 2, y1),
            y2 = ifelse(type2 == "sink", max(y2) / 2, y2),
            y2 = ifelse(type2 == "sink", NA, y2),
            y1 = ifelse(type1 == "growth", y2 + 2.5, y1),
            y2 = ifelse(type2 == "growth", y1 + 2.5, y2),
            y1 = ifelse(type1 == "extend_b" | type1 == "extend_f", y1 + dy_ext, y1),
            y2 = ifelse(type2 == "extend_f" | type2 == "extend_b", y2 + dy_ext, y2),
            mc1 = ifelse(type1 == "src", mc2, mc1)
        ) %>%
        as_tibble()

    nn <- nn %>%
        left_join(metacell_type%>% select(mc1 = metacell, col1 = mc_col, cell_type1 = cell_type), by = c("mc1")) %>%
        left_join(metacell_type%>% select(mc2 = metacell, col2 = mc_col, cell_type2 = cell_type), by = c("mc2"))

    nn <- nn %>%
        mutate(edge_lwd = pmin(flow / edge_w_scale, max_lwd)) %>%
        mutate(edge_type = ifelse(type2 == "growth", "growth", cell_type1))
    # mutate(edge_col = ifelse(type2 == "growth", "black", col1))

    cols <- get_cell_type_colors(dataset)
    cols[["growth"]] <- "black"

    max_time <- length(mc_propagate)
    m1 <- as.numeric(nn$mc1)
    m2 <- as.numeric(nn$mc2)

    head(mc_propagate[[1]]) # This line is there in order to load the required 'Matrix' package

    max_m <- ncol(mc_propagate[[1]])
    prop_flow <- rep(0, nrow(nn))
    for (t in 1:max_time) {
        f <- (nn$time1 == t) & nn$mc1 > 0 & nn$mc2 > 0
        prop_flow[f] <- mc_propagate[[t]][m1[f] + max_m * (m2[f] - 1)]
    }

    prop_flow[prop_flow / edge_w_scale < 1] <- 0

    nn <- nn %>%
        mutate(prop_flow = prop_flow) %>%
        mutate(edge_lwd_fg = pmin(prop_flow / edge_w_scale, max_lwd)) %>%
        filter(!is.na(x1), !is.na(x2), !is.na(y1), !is.na(y2))

    p <- nn %>%
        ggplot(aes(x = x1, y = y1, xend = x2, yend = y2)) +
        geom_segment(data = nn %>% filter(!is.na(edge_lwd)), aes(size = edge_lwd), col = bg_col) +
        geom_segment(data = nn %>% filter(prop_flow > 0, !is.na(edge_type)), aes(col = edge_type, size = edge_lwd_fg)) +
        scale_color_manual(name = "", values = cols) +
        scale_x_continuous(
            limits = c(0, max(nn$time2) + 1),
            labels = function(x) c("", round(time_annot$time_desc, digits = 2), "")[x + 1]
        ) +
        scale_size_continuous(range = c(0, 2)) +
        xlab("Time (E[t])") +
        ylab("") +
        annotate("segment", x = 1, xend = max(nn$time2) - 1, y = as.numeric(mc_rank[as.character(metacell)]), yend = as.numeric(mc_rank[as.character(metacell)]), linetype = "dashed") +
        annotate("text", x = max(nn$time2), y = as.numeric(mc_rank[as.character(metacell)]), label = paste0("MC ", metacell)) +
        theme_classic(15) +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank()) +
        guides(size = "none", color = guide_legend(override.aes = list(size = 5))) +
        theme(legend.position = "bottom")

    return(p)
}
