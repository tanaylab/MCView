load_flow <- function(flow_file, mc_egc, project, dataset, time_annotation_file = NULL, gene_names = NULL) {
    flow <- parse_flow(flow_file, colnames(mc_egc))

    if (!is.null(time_annotation_file)) {
        time_annotation <- parse_time_annotation(time_annotation_file, flow)
    } else {
        time_bins <- unique(c(flow$time1, flow$time2))
        time_annotation <- tibble::tibble(time_bin = time_bins, label = time_bins)
    }

    cache_dir <- project_cache_dir(project)
    serialize_shiny_data(flow, "flow", dataset = dataset, cache_dir = cache_dir)
    serialize_shiny_data(time_annotation, "time_annotation", dataset = dataset, cache_dir = cache_dir)
    cli_alert_info("Loaded flow information. Number of time bins: {.val {nrow(time_annotation)}}")
}

parse_flow <- function(flow, metacell_names) {
    if (is.character(flow)) {
        file <- flow
        flow <- fread(flow) %>% as_tibble()
    } else {
        file <- "flow"
    }

    purrr::walk(c("mc1", "mc2", "time1", "time2", "type1", "type2", "flow"), ~ {
        if (!has_name(flow, .x)) {
            cli_abort("{.field {file}} should have a column named {.field {.x}}")
        }
    })

    missing_mcs <- setdiff(unique(c(flow$mc1, flow$mc2)), metacell_names)
    missing_mcs <- missing_mcs[!(missing_mcs %in% c("src", "sink"))]
    if (length(missing_mcs) > 0) {
        cli_abort("Flow file contains metacells that are missing from the data: {.val {missing_mcs}}")
    }

    missing_mcs_flow <- setdiff(metacell_names, unique(c(flow$mc1, flow$mc2)))
    missing_mcs_flow <- missing_mcs_flow[!(missing_mcs_flow %in% c("src", "sink"))]
    if (length(missing_mcs_flow) > 0) {
        cli_warn("Flow file is missing metacells from the data ({.val {length(missing_mcs_flow)}} metacells): {.val {missing_mcs_flow}}")
    }

    return(flow)
}

parse_time_annotation <- function(time_annotation, flow) {
    if (is.character(time_annotation)) {
        file <- time_annotation
        time_annotation <- fread(time_annotation) %>% as_tibble()
    } else {
        file <- "time_annotation"
    }

    purrr::walk(c("time_bin", "label"), ~ {
        if (!has_name(time_annotation, .x)) {
            cli_abort("{.field {file}} should have a column named {.field {.x}}")
        }
    })
    flow_f <- flow %>%
        filter(mc1 != "src", mc2 != "sink")
    time_bins <- unique(c(flow_f$time1, flow_f$time2))
    missing_time_bins <- setdiff(time_bins, time_annotation$time_bin)
    if (length(missing_time_bins) > 0) {
        cli_abort("Time annotation file contains time bins that are missing from the flow file: {.val {missing_time_bins}}")
    }

    missing_time_bins_flow <- setdiff(time_bins, time_annotation$time_bin)
    if (length(missing_time_bins_flow) > 0) {
        cli_abort("Time annotation file is missing time bins from the flow file: {.val {missing_time_bins_flow}}")
    }

    return(time_annotation)
}
