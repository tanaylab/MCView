#' st_flow UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_st_flow_ui <- function(id) {
    ns <- NS(id)
    tagList(
        generic_column(
            width = 12,
            generic_box(
                title = "Temporal Flow",
                width = 12,
                plotOutput(ns("Temporal_Flow"), height = "500px"),
                status = "primary",
                solidHeader = TRUE
            )
        ),        
        generic_column(
            width = 12,
            generic_box(
                title = "Type composition",
                width = 12,
                plotOutput(ns("Type_composition"), height = "200px"),
                status = "primary",
                solidHeader = TRUE
            )
        )
    )
}


#' st_flow sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_st_flow_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            shinyWidgets::radioGroupButtons(
                    inputId = ns("mode"),
                    label = "Select display:",
                    choices = c(
                        "Types",
                        "SMCs"
                    ),
                    selected = "Types",
                    justified = TRUE
                ),
            uiOutput(ns("display_select")),
            uiOutput(ns("norm_flow"))
        )
    )
}


#' st_flow Server Function
#'
#' @noRd
mod_st_flow_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)
            display_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors)
            norm_selector(input, output, session, dataset, ns)

            data <- get_mc_data(dataset(), "spatial_flow_data")

            flow_to = summarise_flow_to(data$ed)
            flow_from = summarise_flow_from(data$ed)
            flow_to_spat = summarise_flow_to_spatial(data$ed)
            flow_from_spat = summarise_flow_from_spatial(data$ed)

            output$Temporal_Flow = renderPlot({plot_temporal_flow_bars(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors,
                                                                       flow_to, flow_from, flow_to_spat, flow_from_spat)})
            output$Type_composition = renderPlot({plot_type_composition(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors)})
        }
    )
}

plot_type_composition = function(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors){
    
    req(input$mode)

    metacell_types_df = metacell_types()

    if(input$mode == 'Types'){
            req(input$display_select %in% metacell_types_df$cell_type)
            selected = input$display_select
    }else if(input$mode == 'SMCs'){
            req(input$display_select %in% metacell_names())
            smc = input$display_select
            selected = metacell_types_df[metacell_types_df$metacell == smc,]$cell_type
    }

    flow = data$f_sm_sb_tb
    flow = flow[flow$cell_type == selected,] %>% group_by(smc, time_bin) %>% summarise(f = sum(f))
    flow_total = flow %>% group_by(smc) %>% summarise(f = sum(f))
    relevant_smcs = flow_total[flow_total$f > 0,]$smc
    flow = flow[flow$smc %in% relevant_smcs,]
    flow$time_bin = factor(flow$time_bin, levels = as.character(sort(as.numeric(unique(flow$time_bin)))))

    if(input$mode == 'SMCs'){
        flow$selected = flow$smc == smc
    }else{
        flow$selected = FALSE
    }

    ctype_color = metacell_types_df$mc_col
    names(ctype_color) = metacell_types_df$metacell

    g = ggplot(flow, aes(y=f, x=smc, fill = smc)) + 
            geom_bar(stat = "identity", aes(color = selected), linewidth = 1.5) + 
            facet_wrap(~time_bin, nrow = 1) + 
            scale_fill_manual(values=ctype_color) + 
            scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
            guides(color = "none") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                    legend.title = element_blank(),
                    axis.title = element_blank(),   
                    axis.text = element_text(size = 16),
                    legend.text = element_text(size = 18),
                    plot.title = element_text(size = 20, face = "bold"),
                    strip.text = element_text(size = 18)) + ggtitle(selected) + guides(fill="none")

    return(g)

}

get_spat_dict = function(){
        spat_dict = c('M1' = 'Rostral', 'M2' = 'Rostral', 'M3' = 'Rostral', 'M4' = 'Rostral', 
                  'M5' = 'Distal', 'M6' = 'Distal', 'M7' = 'Distal', 'M8' = 'Distal', 
                  'M9' = 'Caudal', 'M10' = 'Caudal', 'M11' = 'Caudal', 'M12' = 'Caudal',
                  'L1' = 'Lateral', 'L2' = 'Lateral', 'L3' = 'Lateral', 'L4' = 'Lateral', 
                  'L5' = 'Lateral', 'L6' = 'Lateral', 'I' = 'Lateral', 
                  'X1_1' = 'X', 'X1_2' = 'X', 'X1_3' = 'X', 'X1_4' = 'X')
        return(spat_dict)
}

summarise_flow_to = function(ed){

    flow_to = ed %>% group_by(time_bin, smc1, smc2) %>% summarize(f = sum(flow))
    flow_to = flow_to %>% group_by(time_bin, smc2) %>% mutate(f_norm = f/sum(f))

    return(flow_to)
}

summarise_flow_to_spatial = function(ed){
    spat_dict = get_spat_dict()

    flow_to = ed
    flow_to$spat1 = spat_dict[flow_to$sbin1]
    flow_to$spat2 = spat_dict[flow_to$sbin2]
    flow_to = flow_to[flow_to$spat1 != 'X',]

    flow_to = flow_to %>% group_by(time_bin, smc1, smc2, spat1, spat2) %>% summarize(f = sum(flow))
    flow_to = flow_to %>% group_by(time_bin, smc2) %>% mutate(f_norm = f/sum(f))

    return(flow_to)
}

summarise_flow_from = function(ed){

    flow_from = ed %>% group_by(time_bin, smc1, smc2) %>% summarize(f = sum(flow))
    flow_from = flow_from %>% group_by(time_bin, smc1) %>% mutate(f_norm = f/sum(f))

    return(flow_from)
}

summarise_flow_from_spatial = function(ed){

    spat_dict = get_spat_dict()

    flow_from = ed
    flow_from$spat1 = spat_dict[flow_from$sbin1]
    flow_from$spat2 = spat_dict[flow_from$sbin2]
    flow_from = flow_from[flow_from$spat2 != 'X',]

    flow_from = flow_from %>% group_by(time_bin, smc1, smc2, spat1, spat2) %>% summarize(f = sum(flow))
    flow_from = flow_from %>% group_by(time_bin, smc1) %>% mutate(f_norm = f/sum(f))

    return(flow_from)
}

plot_temporal_flow_bars = function(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors, flow_to, flow_from, flow_to_spat, flow_from_spat){
    
    req(input$mode)
    req(input$spread_spatial)
    
    if(input$spread_spatial == 'None'){
        flow_to = flow_to
        flow_from = flow_from
    }else{
        flow_to = flow_to_spat
        flow_from = flow_from_spat
    }

    time_bins = sort(as.numeric(unique(flow_to$time_bin)))
    time_bins_sub = time_bins[2:length(time_bins)]
    flow_to = flow_to[flow_to$time_bin %in% as.character(time_bins_sub),]
    flow_from = flow_from[flow_from$time_bin %in% as.character(time_bins_sub),]
    
    metacell_types_df = metacell_types()
    metacell_types_l = metacell_types_df$cell_type
    names(metacell_types_l) = metacell_types_df$metacell

    expand_type = ifelse(is.null(input$Expand_Type), FALSE, input$Expand_Type)

    if(input$mode == 'Types' | (input$mode == 'SMCs' & expand_type)){

        flow_from$ent1 = metacell_types_l[flow_from$smc1]
        flow_from$ent2 = metacell_types_l[flow_from$smc2]

        flow_to$ent1 = metacell_types_l[flow_to$smc1]
        flow_to$ent2 = metacell_types_l[flow_to$smc2]

        if(input$mode == 'Types'){
            req(input$display_select %in% metacell_types_df$cell_type)
            selected = input$display_select
        }else if(input$mode == 'SMCs' & expand_type){
            req(input$display_select %in% metacell_names())
            smc = input$display_select
            selected = metacell_types_df[metacell_types_df$metacell == smc,]$cell_type
        }
        
        if(input$spread_spatial == 'None'){
            flow_from = flow_from %>% group_by(time_bin, ent1, ent2) %>% summarise(f = sum(f, na.rm = T))
            flow_from = flow_from %>% group_by(time_bin, ent1) %>% mutate(f_norm = f/sum(f, na.rm = T))

            flow_to = flow_to %>% group_by(time_bin, ent1, ent2) %>% summarise(f = sum(f, na.rm = T))
            flow_to = flow_to %>% group_by(time_bin, ent2) %>% mutate(f_norm = f/sum(f))
        }else{
            flow_from = flow_from %>% group_by(time_bin, ent1, ent2, spat1, spat2) %>% summarise(f = sum(f, na.rm = T))
            flow_from = flow_from %>% group_by(time_bin, ent1) %>% mutate(f_norm = f/sum(f, na.rm = T))

            flow_to = flow_to %>% group_by(time_bin, ent1, ent2, spat1, spat2) %>% summarise(f = sum(f, na.rm = T))
            flow_to = flow_to %>% group_by(time_bin, ent2) %>% mutate(f_norm = f/sum(f))
        }

        ctype_color = cell_type_colors()$color
        names(ctype_color) = cell_type_colors()$cell_type
        ctype_color[c('source', 'sink')] = 'gray'

    }else if(input$mode == 'SMCs'){
        req(input$display_select %in% metacell_names())
        selected = input$display_select
        flow_from$ent1 = flow_from$smc1
        flow_from$ent2 = flow_from$smc2
        flow_to$ent1 = flow_to$smc1
        flow_to$ent2 = flow_to$smc2

        ctype_color = metacell_types_df$mc_col
        names(ctype_color) = metacell_types_df$metacell
        ctype_color[c('source', 'sink')] = 'gray'
    }

    flow_from[is.na(flow_from$ent2),]$ent2 = 'sink'
    flow_from[is.na(flow_from$ent1),]$ent1 = 'source'
    flow_from[is.na(flow_from$f_norm),]$f_norm = 0

    flow_to[is.na(flow_to$ent2),]$ent2 = 'sink'
    flow_to[is.na(flow_to$ent1),]$ent1 = 'source'
    flow_to[is.na(flow_to$f_norm),]$f_norm = 0
    
    if(input$spread_spatial %in% c("Rostral","Distal","Lateral", "Caudal")){
        subset_flow_from = flow_from[flow_from$ent1 == selected & flow_from$spat1 == input$spread_spatial,]
        subset_flow_to = flow_to[flow_to$ent2 == selected & flow_to$spat2 == input$spread_spatial,]
    }else{
        subset_flow_from = flow_from[flow_from$ent1 == selected,]
        subset_flow_to = flow_to[flow_to$ent2 == selected,]
    }

    subset_flow_from = subset_flow_from[subset_flow_from$f_norm > 0.05 | subset_flow_from$ent2 == selected,]
    subset_flow_from$time_bin = as.character(as.numeric(subset_flow_from$time_bin)-1)

    subset_flow_to = subset_flow_to[subset_flow_to$f_norm > 0.05 | subset_flow_to$ent1 == selected,]

    subset_flow = data.frame(time_bin = c(subset_flow_from$time_bin, subset_flow_to$time_bin),
                             smc = c(subset_flow_from$ent2, subset_flow_to$ent1),
                             f = c(subset_flow_from$f, -subset_flow_to$f),
                             f_norm = c(subset_flow_from$f_norm, -subset_flow_to$f_norm))

    if(!input$spread_spatial == 'None'){
        subset_flow$spat = factor(c(subset_flow_from$spat2, subset_flow_to$spat1), levels = c("Rostral","Distal","Lateral", "Caudal"))
        subset_flow <- subset_flow %>% complete(time_bin, smc, spat = c("Rostral","Distal","Lateral", "Caudal"), fill = list(f = 0, f_norm = 0))
    }

    subset_flow$smc = factor(subset_flow$smc, levels = unique(subset_flow$smc))

    subset_flow$time_bin = factor(subset_flow$time_bin, levels = as.character(time_bins))

    if(input$norm_flow){
       subset_flow$flow_plot = subset_flow$f_norm
    }else if(!input$norm_flow){
       subset_flow$flow_plot = subset_flow$f
    }

    flow_label = paste0(time_bins, "->", time_bins[-1])[-length(time_bins)]
    flow_label = c(paste0('source', "->", time_bins[1]), flow_label, paste0(time_bins[length(time_bins)], "->", 'sink'))

    flow_label = paste0(flow_label[seq(1,(length(flow_label)-1))],  '    ', 
                        flow_label[seq(2,(length(flow_label)))])
    names(flow_label) = as.character(time_bins)

    if(input$mode == 'SMCs'){
        subset_flow$mark = subset_flow$smc == selected
    }else{
        subset_flow$mark = FALSE
    }

    x_lim = max(abs(min(subset_flow$flow_plot, na.rm = T)), max(subset_flow$flow_plot, na.rm = T))

    if(input$spread_spatial %in% c("Rostral","Distal","Lateral", "Caudal")){ # "All", 
        # browser()
        dashes_df = data.frame(smc = as.numeric(unique(subset_flow$smc)))
        sub_dashes = data.frame(spat = c(dashes_df$smc, dashes_df$smc + 0.25, dashes_df$smc + 0.5, dashes_df$smc + 0.75))

        g = ggplot(subset_flow, aes(x = flow_plot, y = smc, fill = smc)) + 
                        geom_bar(stat = "identity", position = 'dodge', aes(color = mark, group = spat), linewidth = 1.5) +
                        facet_wrap(~time_bin, nrow = 1, labeller = labeller(time_bin = flow_label)) + 
                        scale_fill_manual(values=ctype_color) + 
                        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
                        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
                        geom_hline(data = dashes_df, aes(yintercept = smc - 0.5, linetype = "dotted", color = "gray")) +
                        geom_hline(data = sub_dashes, aes(yintercept = spat - 0.5, linetype = "dotted", color = "black")) +
                        xlim(-x_lim, x_lim) + guides(color = "none") +
                        theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                              legend.title = element_blank(),
                              axis.title = element_blank(),   
                              axis.text = element_text(size = 16),
                              legend.text = element_text(size = 18),
                              plot.title = element_text(size = 20, face = "bold"),
                              strip.text = element_text(size = 18)) + ggtitle(selected)
    }else{
        
        g = ggplot(subset_flow, aes(x = flow_plot, y = smc, fill = smc)) + 
                        geom_bar(stat = "identity", aes(color = mark), linewidth = 1.5) +
                        facet_wrap(~time_bin, nrow = 1, labeller = labeller(time_bin = flow_label)) + 
                        scale_fill_manual(values=ctype_color) + 
                        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
                        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
                        xlim(-x_lim, x_lim) + guides(color = "none") +
                        theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                              legend.title = element_blank(),
                              axis.title = element_blank(),   
                              axis.text = element_text(size = 16),
                              legend.text = element_text(size = 18),
                              plot.title = element_text(size = 20, face = "bold"),
                              strip.text = element_text(size = 18)) + ggtitle(selected)
    }
    
    # file_name = selected
    # if(grepl('/', file_name)){
    #     file_name = gsub(x = file_name, pattern = '/', replacement = '_')
    # }
    # if(grepl('\\?', file_name)){
    #     file_name = gsub(x = file_name, pattern = '\\?', replacement = 'o')
    # }
    # ggsave(paste0(file_name, '_', input$spread_spatial, '.png'))

    return(g)
}

display_selectors <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors) {
    output$display_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "SMCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())

            md <- get_mc_data(dataset(), "metadata")
            if (!is.null(md) && has_name(md, "similar")) {
                choices <- metacell_names()
                names(choices) <- ifelse(md$similar == "dissimilar", paste0(metacell_names(), " (dissimilar)"), metacell_names())
            } else {
                choices <- metacell_names()
            }
            tagList(
            shinyWidgets::pickerInput(ns("display_select"), "Smc",
                choices = choices,
                selected = 'M267.82',
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            ),
            shinyWidgets::switchInput(
                inputId = ns("Expand_Type"),
                label = "Show Type",
                size = "sm",
                value = FALSE
            )
            )
        }else if(input$mode == "Types"){
            req(cell_type_colors())
            req(metacell_types())
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            tagList(
                shinyWidgets::pickerInput(ns("display_select"), "Cell Type",
                    choices = cell_types,
                    selected = 'Epiblast',
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
                )
            )
        }
    })}

norm_selector <- function(input, output, session, dataset, ns) {
    output$norm_flow <- renderUI({
        tagList(
         shinyWidgets::pickerInput(
                inputId = ns("spread_spatial"), 
                label = "Spread Spatial",
                choices = c('None', 'Rostral', 'Distal', 'Caudal', 'Lateral'), # 'All', 
                selected = 'None',
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)
            ),
        shinyWidgets::switchInput(
                inputId = ns("norm_flow"),
                label = "Normalize flow",
                size = "sm",
                value = FALSE
            )
        )
        }
)}
