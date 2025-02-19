#' SPATIAL UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_spatial_ui <- function(id) {
    ns <- NS(id)
    tagList(
        generic_column(
            width = 12,
            generic_box(
                title = "A-P probability",
                width = 12,
                plotOutput(ns("APplot"), height = "200px"),
                status = "primary",
                solidHeader = TRUE
            )
        ),
        generic_column(
            width = 12,
            generic_box(
                title = "Spatial Density",
                width = 12,
                plotOutput(ns("beatleplot"), height = "200px"),
                status = "primary",
                solidHeader = TRUE
            )
        )

    )
}

mod_spatial_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            div(
                # id = ns("sidebar_select"),
                shinyWidgets::radioGroupButtons(
                    inputId = ns("mode"),
                    label = "Select display:",
                    choices = c(
                        "Types",
                        "SMCs",
                        "Genes"
                    ),
                    selected = "Types",
                    justified = TRUE
                ),
                uiOutput(ns("display_select")),
                uiOutput(ns("types_express_genes"))
            )
        )
    )
}

mod_spatial_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    library(here)
    library(gridExtra)
    library(grid)
    devtools::load_all(here("mingeom"), export_all = FALSE, quiet = TRUE)
    devtools::load_all(here("spatula"), export_all = FALSE, quiet = TRUE)

    moduleServer(
        id,
        function(input, output, session) {
            ns <- NS(id)            

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)
            display_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors)
            mod_spatial_globals_observers(input, session, globals)

            data <- get_mc_data(dataset(), "spatial_flow_data")

            output$beatleplot <- renderPlot({      
                plot_map(input, data, dataset, metacell_names, metacell_types)
            }, height = function(){200*plot_height(input, data)})

            output$APplot <- renderPlot({
                plot_AP(input, data, dataset, metacell_names, metacell_types)
            })

            output$types_express_genes <- renderUI({
                req(input$mode == 'Genes')
                types_express_genes_selector(input, dataset, ns, "types_exp_genes", metacell_types, cell_type_colors)
            })
 
    })
}

types_express_genes_selector <- function(input, dataset, ns, id, metacell_types, cell_type_colors) {
    req(input$display_select %in% gene_names(dataset()))
    req(metacell_types())
    
    top_types = get_top_egc_types(input, dataset, metacell_types)
    tagList(
        selectInput(
            ns(id),
            label = 'types expressing the gene',
            choices = top_types,
            selected = NULL,
            size = 10,
            selectize = FALSE,
            multiple = TRUE
        )
    )
}

get_top_egc_types <- function(input, dataset, metacell_types) {

    req(input$display_select %in% gene_names(dataset()))
    egc <- get_gene_egc(input$display_select, dataset(), atlas = FALSE) + egc_epsilon
    type_exp = tapply(egc, metacell_types()$cell_type, mean)

    show_types = head(sort(type_exp, decreasing = T), 10)
    show_types = paste0('[', formatC(unname(show_types), format = "e", digits = 1),'] ', names(show_types))
    return(show_types)
}

display_selectors <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors) {
    output$display_select <- renderUI({
        req(dataset())
        req(input$mode)
        if(input$mode == "Types"){
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
                ),
                shinyWidgets::switchInput(
                inputId = ns("Expand_Smcs"),
                label = "Show SMCs",
                size = "sm",
                value = FALSE
                )
            )
        } else if (input$mode == "SMCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            # add 'similar' annotation
            md <- get_mc_data(dataset(), "metadata")
            if (!is.null(md) && has_name(md, "similar")) {
                choices <- metacell_names()
                names(choices) <- ifelse(md$similar == "dissimilar", paste0(metacell_names(), " (dissimilar)"), metacell_names())
            } else {
                choices <- metacell_names()
            }
            shinyWidgets::pickerInput(ns("display_select"), "Smc",
                choices = choices,
                selected = 'M267.82',
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        }  else if (input$mode == "Genes") {
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type

            gene_choices <- gene_names(dataset())
            tagList(
                shinyWidgets::pickerInput(ns("display_select"), "Gene",
                    choices = gene_choices,
                    selected = 'Sfrp1',
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
                ),
                shinyWidgets::switchInput(
                inputId = ns("total_max"),
                label = "global",
                size = "sm",
                value = TRUE
                ),
                shinyWidgets::pickerInput(ns("inspect_sub_manif"), "Inspect Sub Manif",
                    choices = c(c('ALL', 'Endoderm', 'Mesoderm', 'Ectoderm'), cell_types),
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", c(c('#000000', '#e68ce3', '#28aaeb', '#3dd11b'), cell_types_hex), ";")
                    )
                )
            )
        }
    })}

mod_spatial_globals_observers <- function(input, session, globals, notification_suffix = " in \"Spatial\" tab") {
    observe({
        req(input$mode == "SMCs")
        req(globals$selected_smc)
        req(input$display_select)
        shinyWidgets::updatePickerInput(session, "display_select", selected = globals$selected_smc)
        showNotification(glue("Selected {globals$selected_smc}{notification_suffix}"))
        globals$selected_smc <- NULL
    })
}

plot_height <- function(input, data){
    ctype_manif = c()
    if(input$mode == 'Types' & input$Expand_Smcs){
       to_plot = input$display_select
       ctype_manif = unique(data$f_sm_sb_tb$smc[data$f_sm_sb_tb$cell_type == to_plot])
    }
    return(length(ctype_manif) + 1)
}

summarise_flow <- function(input, data, mode, tbin_time){

        to_plot = input$display_select

        smc = NULL
        ctype_manif = NULL
        # P(s|t)
        f_sb_tb = data$f_sm_sb_tb %>% group_by(sbin, time_bin) %>% summarise(p_s = sum(f))

        if(mode == 'Types'){
            # P(ctype|t)
            f_ct_sb_tb = data$f_sm_sb_tb %>% group_by(cell_type, sbin, time_bin, n_col, layer, age, p_s_t_tbin_model_input) %>% 
                                            summarise(f_ctype = sum(f), n_th = sum(n_th))

            f_ct_sb_tb = f_ct_sb_tb %>% group_by(cell_type, time_bin) %>% 
                                    mutate(p_ct_sb = f_ctype/sum(f_ctype))
        }else if(mode == 'SMCs'){
            smc = T
            f_ct_sb_tb = data$f_sm_sb_tb

            f_ct_sb_tb = f_ct_sb_tb %>% group_by(smc, time_bin) %>% 
                                    mutate(p_ct_sb = f/sum(f))
        }

        f_ct_sb_tb = left_join(f_ct_sb_tb, f_sb_tb, by = c('sbin', 'time_bin'))

        f_ct_sb_tb$density = f_ct_sb_tb$p_ct_sb/f_ct_sb_tb$p_s

        return(f_ct_sb_tb)
}

summarise_insitu = function(input, data, dataset, mode = 'egc'){

    f_ct_sb_tb = data$f_sm_sb_tb
    sub_manif = unique(f_ct_sb_tb$cell_type)

    layer = input$inspect_sub_manif
    if(layer == 'ALL'){
        sub_manif = unique(f_ct_sb_tb$cell_type)
    }else if(layer == 'Endoderm'){
        sub_manif = unique(f_ct_sb_tb[f_ct_sb_tb$layer == 'Endo', ]$cell_type)
    }else if(layer == 'Mesoderm'){
        sub_manif = unique(f_ct_sb_tb[f_ct_sb_tb$layer == 'Meso', ]$cell_type)
    }else if(layer == 'Ectoderm'){
        sub_manif = unique(f_ct_sb_tb[f_ct_sb_tb$layer == 'Ecto', ]$cell_type)
    }else{
        sub_manif = unique(f_ct_sb_tb[f_ct_sb_tb$cell_type == layer, ]$cell_type)
    }

    f_ct_sb_tb = f_ct_sb_tb[f_ct_sb_tb$cell_type %in% sub_manif, ]
    # f_ct_sb_tb = f_ct_sb_tb %>% group_by(time_bin) %>% mutate(layer_frac = sum(f))
    f_ct_sb_tb = f_ct_sb_tb %>% group_by(time_bin, sbin, smc) %>% mutate(f_m_s = sum(f))
    f_ct_sb_tb = f_ct_sb_tb %>% group_by(time_bin, sbin) %>% mutate(f_s = sum(f))

    egc <- get_gene_egc(input$display_select, dataset(), atlas = FALSE) + egc_epsilon
    f_ct_sb_tb$egc = egc[as.character(f_ct_sb_tb$smc)]
    tot_egc = f_ct_sb_tb %>% group_by(time_bin) %>% summarise(tot_egc = sum(f_m_s*egc, na.rm = T))

    if(mode == 'dens'){
        spat_density = f_ct_sb_tb %>% group_by(time_bin, sbin) %>% summarise(density = sum(f_m_s*egc/f_s, na.rm = T))
        return(list(spat_density = spat_density, tot_egc = tot_egc))
    }else if(mode == 'egc'){
        spat_egc = f_ct_sb_tb %>% group_by(time_bin, sbin) %>% summarise(spat_egc = sum(f_m_s*egc))
        return(list(spat_egc = spat_egc, tot_egc = tot_egc))
    }
}


plot_AP <- function(input, data, dataset, metacell_names, metacell_types){

    if(input$mode == 'Types'){
        req(input$display_select %in% metacell_types()$cell_type)
    }else if(input$mode == 'SMCs'){
        req(input$display_select %in% metacell_names())
    }else if(input$mode == 'Genes'){
        req(input$display_select %in% gene_names(dataset()))
    }

    to_plot = input$display_select

    if(input$mode == 'Types' | input$mode == 'SMCs'){
        f_ct_sb_tb = summarise_flow(input, data, mode = input$mode, tbin_time)

        if(input$mode == 'Types'){
                plot_data = f_ct_sb_tb[f_ct_sb_tb$cell_type == to_plot & f_ct_sb_tb$sbin %in% paste0('M', seq(1,12)),]
        }else if(input$mode == 'SMCs'){
                plot_data = f_ct_sb_tb[f_ct_sb_tb$smc == to_plot & f_ct_sb_tb$sbin %in% paste0('M', seq(1,12)),]
        }

        plot_data$sbin = factor(plot_data$sbin, levels = paste0('M', seq(1,12)))
        plot_data$time_bin = factor(plot_data$time_bin, levels = as.character(sort(as.integer(unique(plot_data$time_bin)))))

        plot_data$color_plot = 'enough data'
        plot_data[plot_data$n_col < data$th_col,'color_plot'] = 'too little colored cells'
        plot_data[plot_data$n_th < data$th_n,'color_plot'] = 'too little cells'

        lable_colors = c("enough data" = "black", 'too little colored cells'="#969595", 'too little cells'="#d6d4d4")

        p = ggplot(data = plot_data, 
                    aes(x = sbin, y = p_ct_sb, group=1, color = color_plot, size = 1)) + geom_line() +
                    scale_color_manual(values=lable_colors) + ylim(c(egc_epsilon,NA))+
                    facet_wrap(~time_bin, nrow = 1) + ylab('') + 
                    ggtitle(sprintf('P(sbin | %s)', to_plot)) + theme(legend.position="none")
        
        return(p)

    }else if(input$mode == 'Genes'){

        summarise_insitu_ret = summarise_insitu(input, data, dataset)
        spat_egc = summarise_insitu_ret$spat_egc
        tot_egc = summarise_insitu_ret$tot_egc

        plot_data = spat_egc[spat_egc$sbin %in% paste0('M', seq(1,12)),]
        plot_data$sbin = factor(plot_data$sbin, levels = paste0('M', seq(1,12)))
        plot_data$time_bin = factor(plot_data$time_bin, levels = as.character(sort(as.integer(unique(plot_data$time_bin)))))

        p = ggplot(data = plot_data, 
                    aes(x = sbin, y = spat_egc, group=1, size = 1)) + geom_line() +
                    facet_wrap(~time_bin, nrow = 1) + ylab('') + 
                    ggtitle(sprintf('P(sbin | %s)', to_plot)) + theme(legend.position="none")
        # plot_data = spat_egc[spat_egc$cell_type == to_plot & f_ct_sb_tb$sbin %in% paste0('M', seq(1,12)),]

        return(p)
    }
}

plot_map <- function(input, data, dataset, metacell_names, metacell_types) {
    
    if(input$mode == 'Types'){
        req(input$display_select %in% metacell_types()$cell_type)
    }else if(input$mode == 'SMCs'){
        req(input$display_select %in% metacell_names())
    }else if(input$mode == 'Genes'){
        req(input$display_select %in% gene_names(dataset()))
    }

    smc_idx = 0
    g = NULL

    tbin_time = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$age
    names(tbin_time) = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$time_bin
    tbs = sort(as.integer(unique(data$f_sm_sb_tb$time_bin)))
    to_plot = input$display_select

    p = list()
    tb_idx = 1

    if(input$mode == 'Types' | input$mode == 'SMCs'){

        f_ct_sb_tb = summarise_flow(input, data, mode = input$mode, tbin_time)
        
        for(tb in tbs){

            if(input$mode == 'Types'){
                plot_data = f_ct_sb_tb[f_ct_sb_tb$cell_type == to_plot & f_ct_sb_tb$time_bin == tb,]
            }else if(input$mode == 'SMCs'){
                plot_data = f_ct_sb_tb[f_ct_sb_tb$smc == to_plot & f_ct_sb_tb$time_bin == tb,]
            }
            p[[tb_idx]] = beatle_plot_wrap(input, data, plot_data, tbin_time, tb, smc = input$mode == 'SMCs', f_ct_sb_tb)
            tb_idx = tb_idx + 1
        }

        if(input$mode == 'Types' & input$Expand_Smcs){
            smc_idx = 1
            smc = T
            ctype_manif = unique(data$f_sm_sb_tb$smc[data$f_sm_sb_tb$cell_type == to_plot])

            f_ct_sb_tb = summarise_flow(input, data, mode = 'SMCs', tbin_time)

            smc_idx = 1
            for(smc in ctype_manif){
                tb_idx = 1
                for(tb in tbs){

                    plot_data = f_ct_sb_tb[f_ct_sb_tb$smc == smc & f_ct_sb_tb$time_bin == tb,]
                    p[[length(tbs)*smc_idx + tb_idx]] = beatle_plot_wrap(input, data, plot_data, tbin_time, tb, smc, f_ct_sb_tb, sub_title = T)
                    tb_idx = tb_idx + 1
                }
                smc_idx = smc_idx + 1
            }
        }

    }else if(input$mode == 'Genes'){
    
        summarise_insitu_res= summarise_insitu(input, data, dataset, mode = 'dens')
        spat_density = summarise_insitu_res$spat_density
        egc = summarise_insitu_res$tot_egc

        if(input$total_max == TRUE){
            max_col = max(spat_density$density)
        }else{
            max_col = FALSE
        }

        for(tb in tbs){
            plot_data = spat_density[spat_density$time_bin == tb,]
            plot_data$egc = formatC(egc[egc$time_bin == tb,]$tot_egc, format = "e", digits = 2)

            p[[tb_idx]] = beatle_plot_wrap(input, data, plot_data, tbin_time, tb, smc = input$mode == 'SMCs', f_ct_sb_tb, insitu = T, max_col = max_col)
            tb_idx = tb_idx + 1
        }

    }

    g = arrangeGrob(grobs = lapply(p, ggplotGrob), ncol = length(tbs), top = to_plot)
    return(grid.draw(g))
}

gg_theme = function(){
    return(theme(axis.ticks = element_blank(),
                 axis.text = element_blank(),
                 axis.title=element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank()))
}

beatle_plot_wrap = function(input, data, plot_data, tbin_time, tb, smc, f_ct_sb_tb, sub_title = F, insitu = F, max_col = F){
    if(nrow(plot_data) == 0){
        plot_data = data.frame(sbin = unique(f_ct_sb_tb$sbin))
        plot_data$density = 0
        plot_data$n_th = 0
        plot_data$n_col = 0
    }
    plot_data$density[is.na(plot_data$density)] = 0

    if(insitu){
        color_sacle = c("white", "yellow", "orange", "red")
        title = paste('tb: ', tb, '\n egc: ', format(unique(plot_data$egc), digits=2))
        min_scale = egc_epsilon
    }else{
        min_scale = 0
        if(unique(plot_data$n_th) < data$th_n){
            color_sacle = c("white", "white")
        }else if(unique(plot_data$n_col) < data$th_col){
            color_sacle = c("white", "gray")
        }else{
            if(smc != F){
                color_sacle = c("white","red")
            }else{
                color_sacle = c("white","blue")
            }
        }

        if(sub_title){
            title = paste0('smc:', unique(plot_data$smc), ' [', unique(plot_data$n_th),'n]')
        }else{
            title = paste(tb, '~', tbin_time[as.character(tb)],  
                                    '[',unique(plot_data$n_col), 'c /', unique(plot_data$n_th),'n]')
        }
    }

    if(max_col){
        max_scale = max_col
    }else{
        max_scale = max(plot_data$density)
    }
    p = beatle_plot(data, 
                    to_plot = plot_data[,c('sbin', 'density')], 
                    tb,
                    min_scale = min_scale,
                    max_scale = max_scale,
                    color_scale = color_sacle,
                    title = title)
    return(p)
}

beatle_plot = function(data, to_plot, tb, min_scale, max_scale, color_scale, title){

        spat_model = round(data$tbin_time[[as.character(tb)]],1)
        proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

        colnames(to_plot) = c('sbin', 'f')

        p = proj$heatmap(to_plot, show_poles = TRUE, color = "gray75", size = 0.25) +
                         scale_fill_gradientn(colours = color_scale, name = 'density',
                         limits=c(min_scale, max_scale)) + ggtitle(title) + gg_theme()

        return(p)
}
