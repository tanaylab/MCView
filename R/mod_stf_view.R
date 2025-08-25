#' STF VIEW UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_Stf_view_ui <- function(id) {
    ns <- NS(id)
    tagList(
        generic_column(
            width = 12,
                generic_box(
                    title = "Type composition",
                    width = 12,
                    plotOutput(ns("Type_composition_Stf_view"), height = "200px"),
                    status = "primary",
                    solidHeader = TRUE
                )
        ),
        generic_column(
            width = 12,
            generic_box(
                title = "spatio-temporal flow, spatial view",
                width = 12,
                plotOutput(ns("Stf_view_focus"), height = "230px"),
                status = "primary",
                solidHeader = TRUE
            )
        ),
        generic_column(
            width = 12,
            generic_box(
                width = 12,
                title = "spatio-temporal contributions, spatial view",
                plotOutput(ns("Stf_view_contribs"), height = "300px"),
                status = "primary",
                solidHeader = TRUE
            )
        )
    )
}


mod_Stf_view_sidebar_ui <- function(id) {
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
                        "SMCs"
                    ),
                    selected = "Types",
                    justified = TRUE
                ),
                shinyWidgets::radioGroupButtons(
                    inputId = ns("from_to"),
                    label = "Select direction:",
                    choices = c(
                        "Income",
                        "Outgo"
                    ),
                    selected = "Income",
                    justified = TRUE
                ),
                uiOutput(ns("display_select_type_smc")),
                uiOutput(ns("display_select_sbin"))
            )
        )
    )
}


mod_Stf_view_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    library(here)
    library(gridExtra)
    library(grid)
    devtools::load_all(here("mingeom"), export_all = FALSE, quiet = TRUE)
    devtools::load_all(here("spatula"), export_all = FALSE, quiet = TRUE)
    library(sf)

    moduleServer(
        id,
        function(input, output, session) {
            ns <- NS(id)            

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)
            display_selector_spat_type(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors)
            data <- get_mc_data(dataset(), "spatial_flow_data")

            eps = 1e-6
            data$ed = data$ed[data$ed$flow > eps,]
            
            sbins = unique(data$f_sm_sb_tb$sbin)
            display_selector_sbin(input, output, session, sbins, ns)

            output$Type_composition_Stf_view = renderPlot({plot_type_composition_beatle_flow(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors)})

            output$Stf_view_focus <- renderPlot({      
                g = strip_beatle_flow_plot(input, data, dataset, metacell_names, metacell_types, cell_type_colors, main = T)
                grid.draw(g)
            })

            f_th = 0.1
            output$Stf_view_contribs <- renderPlot({     
                g = strip_beatle_flow_plot(input, data, dataset, metacell_names, metacell_types, cell_type_colors, f_th, main = F)
                grid.draw(g)
            }, height = function(){250*n_contribs_tbs(input, data, metacell_names, metacell_types, f_th)})

    })
}


n_contribs_tbs = function(input, data, metacell_names, metacell_types, f_th){
    
    metacell_types_l = metacell_types()$cell_type
    names(metacell_types_l) = metacell_types()$metacell
    data$ed$type1 = metacell_types_l[data$ed$smc1]
    data$ed$type2 = metacell_types_l[data$ed$smc2]
    
    if(input$mode == 'Types'){
        req(input$display_select_type_smc %in% metacell_types()$cell_type)
    }else if(input$mode == 'SMCs'){
        req(input$display_select_type_smc %in% metacell_names())
    }

    tbin_time = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$age
    names(tbin_time) = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$time_bin

    to_plot = input$display_select_type_smc

    sb = input$display_select_sbin

    flow_plot = summarise_flow_spatial_tbs(data, to_plot, sb, metacell_types, metacell_names, f_th, from = input$from_to == 'Outgo')

    if(input$from_to == 'Outgo'){
        types = unique(flow_plot$plot_dynamics$ent2)
        return(length(types[!is.na(types)]))
    }else{
        types = unique(flow_plot$plot_dynamics$ent1)
        return(length(types[!is.na(types)]))
    }
    
}

strip_beatle_flow_plot = function(input, data, dataset, metacell_names, metacell_types, cell_type_colors, f_th = 0.01, main = F){
    
    metacell_types_l = metacell_types()$cell_type
    names(metacell_types_l) = metacell_types()$metacell
    data$ed$type1 = metacell_types_l[data$ed$smc1]
    data$ed$type2 = metacell_types_l[data$ed$smc2]
    
    if(input$mode == 'Types'){
        req(input$display_select_type_smc %in% metacell_types()$cell_type)
    }else if(input$mode == 'SMCs'){
        req(input$display_select_type_smc %in% metacell_names())
    }
    g = NULL

    tbin_time = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$age
    names(tbin_time) = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$time_bin

    to_plot = input$display_select_type_smc

    if(input$mode == 'Types'){
        cell_type_colors_tab = cell_type_colors()
        color = cell_type_colors_tab[cell_type_colors_tab$cell_type == to_plot,]$color
    }else if(input$mode == 'SMCs'){
        metacell_types_tab = metacell_types()
        color = metacell_types_tab[metacell_types_tab$metacell == to_plot,]$mc_col
    }

    sb = input$display_select_sbin

    flow_plot = summarise_flow_spatial_tbs(data, to_plot, sb, metacell_types, metacell_names, f_th, from = input$from_to == 'Outgo')
    plot_focus = flow_plot$plot_focus[,c('sbin2', 'flow', 'time_bin')]
    plot_dynamics = flow_plot$plot_dynamics
    tbs = as.character(sort(as.numeric(unique(data$ed$time_bin))))

    if(main){
        

        if(input$from_to == 'Outgo'){
            tbs_dynamics = tbs[1:(length(tbs)-1)]
        }else{
            tbs_dynamics = tbs[2:length(tbs)]
        }

        p = list()
        idx = 1

        tbin_time = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$age
        names(tbin_time) = unique(data$f_sm_sb_tb[,c('time_bin', 'age')])$time_bin
        tbin_time_l = paste0(names(tbin_time), ' ~E',tbin_time)
        names(tbin_time_l) = names(tbin_time)

        for(tb in tbs){
        
            bp = scatter_beatle_plot(input, data, tbin_time, tb, 
                                     ungroup(plot_focus[plot_focus$time_bin == tb, c('sbin2', 'flow')]), 
                                     color)

            # if(tb %in% tbs_dynamics){
            #     if(input$from_to == 'Outgo'){
            #         bp = flow_arrow_plot(bp, data, tb, 
            #                          ungroup(plot_dynamics[plot_dynamics$time_bin == as.character(as.numeric(tb)+1), c('ent1', 'ent2', 'sbin1', 'sbin2', 'flow')]))  
            #     }else{
            #         bp = flow_arrow_plot(bp, data, tb, 
            #                          ungroup(plot_dynamics[plot_dynamics$time_bin == tb, c('ent1', 'ent2', 'sbin1', 'sbin2', 'flow')]))     
            #     }
                
            # }

            p[[idx]] = bp + gg_theme() + ggtitle(sprintf('time bin %s ~E%s %s', tb, unname(tbin_time_l[tb]), ifelse(sb == 'ALL', '', sprintf('\n %s', sb))))
            idx = idx + 1
        }
        g = arrangeGrob(grobs = lapply(p, ggplotGrob), nrow = 1)
        # ggsave(file="whatever.pdf", g) #saves g

        return(g)


    }else{

        p = list()
        tb_idx = 1
        type_idx = 0

        if(input$from_to == 'Outgo'){
            to_flows = unique(plot_dynamics$ent2)
            to_flows = c(to_plot, to_flows[to_flows != to_plot])
            to_flows = to_flows[!is.na(to_flows)]

            tbs_dynamics = tbs[1:(length(tbs)-1)]

            for(tb in tbs_dynamics){
                spat_model = round(data$tbin_time[[as.character(tb)]],1)
                proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

                flow_frac = sum(plot_focus[plot_focus$time_bin == tb,]$flow)

                for(to_flow in to_flows){
                    # "to flow" related to focus
                    plot_focus_to_flow = plot_dynamics[plot_dynamics$ent2 == to_flow & plot_dynamics$time_bin == as.character(as.numeric(tb)+1),]
                    
                    if(sum(plot_focus_to_flow$flow)/flow_frac >= f_th & flow_frac >= 1e-4){

                        if(input$mode == 'Types'){
                            cell_type_colors_tab = cell_type_colors()
                            color2 = cell_type_colors_tab[cell_type_colors_tab$cell_type == to_flow,]$color
                        }else if(input$mode == 'SMCs'){
                            metacell_types_tab = metacell_types()
                            color2 = metacell_types_tab[metacell_types_tab$metacell == to_flow,]$mc_col
                        }

                        plot_focus_from_scatter = plot_focus_to_flow %>% group_by(sbin1) %>% summarise(flow = sum(flow)) 

                        bp = scatter_beatle_plot(input, data, tbin_time, tb, ungroup(plot_focus_from_scatter), color2)
                        if(tb %in% tbs_dynamics){
                            bp = flow_arrow_plot(bp, data, tb, plot_focus_to_flow[,c('ent1', 'ent2', 'sbin1', 'sbin2', 'flow')])
                        }
                        bp = bp + gg_theme() +
                            ggtitle(sprintf('[%s]->\n%s \n time bin %s -> %s %s', 
                                        to_plot, to_flow, tb, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                        p[[tb_idx + type_idx*length(tbs)]] = bp

                    }else{

                        p[[tb_idx + type_idx*length(tbs)]] = proj$heatmap(fill = "white") + gg_theme() +
                            ggtitle(sprintf('[%s]->\n%s \n time bin %s -> %s %s', 
                                        to_plot, to_flow, tb, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    }
                    type_idx = type_idx + 1
                }
                type_idx = 0
                tb_idx = tb_idx + 1
            }

            tb = tbs[length(tbs)]
            spat_model = round(data$tbin_time[[as.character(tb)]],1)
            proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

            for(to_flow in to_flows){
                p[[tb_idx + type_idx*length(tbs)]] = proj$heatmap(fill = "white") + gg_theme() +
                            ggtitle(sprintf('[%s]->\n%s \n time bin %s -> %s %s', 
                                        to_plot, to_flow, tb, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 

                type_idx = type_idx + 1

            }


            g = arrangeGrob(grobs = lapply(p, ggplotGrob), nrow = length(to_flows),
                            top = sprintf('%s %s', to_plot, ifelse(sb == 'ALL', '', sb)))

            return(g)
            
        }else{
            from_flows = unique(plot_dynamics$ent1)
            from_flows = c(to_plot, from_flows[from_flows != to_plot])
            from_flows = from_flows[!is.na(from_flows)]

            tb = tbs[1]
            spat_model = round(data$tbin_time[[as.character(tb)]],1)
            proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

            for(from_flow in from_flows){
                
                p[[tb_idx + type_idx*length(tbs)]] = proj$heatmap(fill = "white") + gg_theme() +
                            ggtitle(sprintf('%s->\n[%s] \n time bin %s -> %s %s', 
                                        from_flow, to_plot, as.character(as.numeric(tb)-1), tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                type_idx = type_idx + 1

            }
            type_idx = 0
            tb_idx = tb_idx + 1
            
            tbs_dynamics = tbs[2:length(tbs)]
                
            for(tb in tbs_dynamics){

                spat_model = round(data$tbin_time[[as.character(tb)]],1)
                proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

                flow_frac = sum(plot_focus[plot_focus$time_bin == tb,]$flow)

                for(from_flow in from_flows){
                    # "to flow" related to focus
                    plot_focus_from_flow = plot_dynamics[plot_dynamics$ent1 == from_flow & plot_dynamics$time_bin == tb,]

                    if(sum(plot_focus_from_flow$flow)/flow_frac >= f_th & flow_frac >= 1e-4){

                        if(input$mode == 'Types'){
                            cell_type_colors_tab = cell_type_colors()
                            color2 = cell_type_colors_tab[cell_type_colors_tab$cell_type == from_flow,]$color
                        }else if(input$mode == 'SMCs'){
                            metacell_types_tab = metacell_types()
                            color2 = metacell_types_tab[metacell_types_tab$metacell == from_flow,]$mc_col
                        }

                        plot_focus_from_scatter = plot_focus_from_flow %>% group_by(sbin1) %>% summarise(flow = sum(flow)) 

                        bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus_from_scatter, color2)
                        if(tb %in% tbs_dynamics){
                            bp = flow_arrow_plot(bp, data, tb, plot_focus_from_flow[,c('ent1', 'ent2', 'sbin1', 'sbin2', 'flow')])
                        }
                        bp = bp + gg_theme() +
                            ggtitle(sprintf('%s->\n[%s] \n time bin %s -> %s %s', 
                                        from_flow, to_plot, as.character(as.numeric(tb)-1), tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                        p[[tb_idx + type_idx*length(tbs)]] = bp
                    }else{
                        p[[tb_idx + type_idx*length(tbs)]] = proj$heatmap(fill = "white") + gg_theme() +
                            ggtitle(sprintf('%s->\n[%s] \n time bin %s -> %s %s', 
                                        from_flow, to_plot, as.character(as.numeric(tb)-1), tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    }
                    type_idx = type_idx + 1
                }
                type_idx = 0
                tb_idx = tb_idx + 1
            }

            g = arrangeGrob(grobs = lapply(p, ggplotGrob), nrow = length(from_flows),
                            top = sprintf('%s %s %s', to_plot, tb, ifelse(sb == 'ALL', '', sb)))

            return(g)
        }
    }
}

summarise_flow_spatial_tbs = function(data, to_plot, sb, metacell_types, metacell_names, f_th, from = TRUE){

    flow = data$ed

    if(to_plot %in% metacell_types()$cell_type){

        flow = flow[,c('type1', 'type2', 'sbin1', 'sbin2', 'time_bin','flow')]
        colnames(flow) = c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin','flow')

    }else if(to_plot %in% metacell_names()){

        flow = flow[,c('smc1', 'smc2', 'sbin1', 'sbin2', 'time_bin','flow')]
        colnames(flow) = c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin','flow')

    }

    if(length(sb) == 1){
        if(sb == 'ALL'){
            sb = unique(flow$sbin2)
        }
    }

    min_th = 0.0001
    plot_focus = flow[flow$ent2 == to_plot & flow$sbin2 %in% sb, ]
    plot_focus = plot_focus %>% group_by(ent2, sbin2, time_bin) %>% summarise(flow = sum(flow))
    f = plot_focus %>% group_by(time_bin) %>% summarise(flow_tot = sum(flow))
    plot_focus = left_join(plot_focus, f, by = c('time_bin'))
    plot_focus = plot_focus[plot_focus$flow_tot > min_th, c('ent2', 'sbin2', 'time_bin', 'flow')]

    if(from){

        plot_dynamics = flow[flow$ent1 == to_plot & flow$sbin1 %in% sb, ]
        plot_dynamics = plot_dynamics %>% group_by(ent1, ent2, sbin1, sbin2, time_bin) %>% summarise(flow = sum(flow))

        f = plot_dynamics %>% group_by(time_bin) %>% summarise(flow_tot = sum(flow))
        plot_dynamics = left_join(plot_dynamics, f, by = c('time_bin'))
        plot_dynamics = plot_dynamics[plot_dynamics$flow_tot > min_th, c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin', 'flow')]

        f = plot_dynamics %>% group_by(time_bin, ent1, ent2) %>% summarise(flow_frac = sum(flow)) %>% 
            group_by(time_bin, ent1) %>% mutate(flow_frac = flow_frac/sum(flow_frac))
        plot_dynamics = left_join(plot_dynamics, f, by = c('time_bin', 'ent1', 'ent2'))
        plot_dynamics = plot_dynamics[plot_dynamics$flow_frac > f_th, c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin', 'flow')]

    }else{

        plot_dynamics = flow[flow$ent2 == to_plot & flow$sbin2 %in% sb, ]
        plot_dynamics = plot_dynamics %>% group_by(ent1, ent2, sbin1, sbin2, time_bin) %>% summarise(flow = sum(flow))

        f = plot_dynamics %>% group_by(time_bin) %>% summarise(flow_tot = sum(flow))
        plot_dynamics = left_join(plot_dynamics, f, by = c('time_bin'))
        plot_dynamics = plot_dynamics[plot_dynamics$flow_tot > min_th, c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin', 'flow')]

        f = plot_dynamics %>% group_by(time_bin, ent1, ent2) %>% summarise(flow_frac = sum(flow)) %>% 
            group_by(time_bin, ent2) %>% mutate(flow_frac = flow_frac/sum(flow_frac))
        plot_dynamics = left_join(plot_dynamics, f, by = c('time_bin', 'ent1', 'ent2'))
        plot_dynamics = plot_dynamics[plot_dynamics$flow_frac > f_th, c('ent1', 'ent2', 'sbin1', 'sbin2', 'time_bin', 'flow')]

    }

    return(list(plot_focus = plot_focus, plot_dynamics = plot_dynamics))
}

display_selector_spat_type <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors) {
    output$display_select_type_smc <- renderUI({
        req(dataset())
        req(input$mode)
        if(input$mode == "Types"){
            req(cell_type_colors())
            req(metacell_types())
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            tagList(
                shinyWidgets::pickerInput(ns("display_select_type_smc"), "Cell Type",
                    choices = cell_types,
                    selected = 'Epiblast',
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
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
            shinyWidgets::pickerInput(ns("display_select_type_smc"), "Smc",
                choices = choices,
                selected = 'M267.82',
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        }
    })}


    display_selector_sbin <- function(input, output, session, sbins, ns) {
        output$display_select_sbin <- renderUI({
            tagList(
                shinyWidgets::pickerInput(ns("display_select_sbin"), "spatial bin",
                    selected = 'ALL',
                    choices = c('ALL', sbins),
                    multiple = TRUE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                )
            )
    })}