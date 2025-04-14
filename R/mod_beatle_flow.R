#' BEATLE FLOW UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_beatle_flow_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 5,
                generic_box(
                    title = "spatio-temporal flow, spatial view",
                    plotOutput(ns("st_flow_spat"), height = "300px"),
                    status = "primary",
                    solidHeader = TRUE
                )
            ),
            generic_column(
                width = 12,
                generic_box(
                    title = "spatio-temporal contributions, spatial view",
                    plotOutput(ns("st_flow_spat_contribs"), height = "300px"),
                    status = "primary",
                    solidHeader = TRUE
                )
            )
        )
    )
}


mod_beatle_flow_sidebar_ui <- function(id) {
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
                        "From",
                        "To"
                    ),
                    selected = "To",
                    justified = TRUE
                ),
                uiOutput(ns("display_select_type_smc")),
                uiOutput(ns("display_select_tbin")),
                uiOutput(ns("display_select_sbin"))
            )
        )
    )
}


mod_beatle_flow_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
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

            eps = 1e-7
            data$ed = data$ed[data$ed$flow > eps,]
            
            tbs = sort(as.numeric(unique(data$f_sm_sb_tb$time_bin)))
            display_selector_tbin(input, output, session, tbs, ns)

            sbins = unique(data$f_sm_sb_tb$sbin)
            display_selector_sbin(input, output, session, sbins, ns)

            output$st_flow_spat <- renderPlot({      
                beatle_flow_plot(input, data, dataset, metacell_names, metacell_types, cell_type_colors, main = T)
            })

            f_th = 0.1
            # browser()
            output$st_flow_spat_contribs <- renderPlot({      
                beatle_flow_plot(input, data, dataset, metacell_names, metacell_types, cell_type_colors, f_th)
            }, height = function(){300*n_contribs(input, data, metacell_names, metacell_types, f_th)})

    })
}

summarise_flow_spatial = function(data, to_plot, tb, sb, metacell_types, metacell_names, from = TRUE){

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

    plot_focus = flow[flow$ent2 == to_plot & flow$sbin2 %in% sb & flow$time_bin == tb, ]
    plot_focus = plot_focus %>% group_by(ent2, sbin2) %>% summarise(flow = sum(flow))

    if(from){

        plot_dynamics = flow[flow$ent1 == to_plot & flow$sbin1 %in% sb & flow$time_bin == as.character(as.numeric(tb)+1), ]
        plot_dynamics = plot_dynamics %>% group_by(ent1, ent2, sbin1, sbin2) %>% summarise(flow = sum(flow))
        
    }else{

        plot_dynamics = flow[flow$ent2 == to_plot & flow$sbin2 %in% sb & flow$time_bin == tb, ]
        plot_dynamics = plot_dynamics %>% group_by(ent1, ent2, sbin1, sbin2) %>% summarise(flow = sum(flow))

    }

    return(list(plot_focus = plot_focus, plot_dynamics = plot_dynamics))
}

n_contribs = function(input, data, metacell_names, metacell_types, f_th){
    
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

    tb = input$display_select_tbin
    sb = input$display_select_sbin

    flow_plot = summarise_flow_spatial(data, to_plot, tb, sb, metacell_types, metacell_names, from = input$from_to == 'From')

    if(input$from_to == 'From'){
        to_flows = unique(flow_plot$plot_dynamics$ent2)

        idx = 0
        for(to_flow in to_flows){
            plot_focus_to_flow = flow_plot$plot_dynamics[flow_plot$plot_dynamics$ent2 == to_flow,] %>% group_by(sbin2) %>% summarise(flow = sum(flow))
            if(sum(plot_focus_to_flow$flow)/sum(flow_plot$plot_focus$flow) >= f_th){
                idx = idx + 1
            }
        }
    }else{

        from_flows = unique(flow_plot$plot_dynamics$ent1)

        idx = 0
        for(from_flow in from_flows){
            plot_focus_from_flow = flow_plot$plot_dynamics[flow_plot$plot_dynamics$ent1 == from_flow,] %>% group_by(sbin1) %>% summarise(flow = sum(flow))
            if(sum(plot_focus_from_flow$flow)/sum(flow_plot$plot_focus$flow) >= f_th){
                idx = idx + 1
            }
        }
    }
    # browser()

    return(idx)
}

beatle_flow_plot = function(input, data, dataset, metacell_names, metacell_types, cell_type_colors, f_th, main = F){
    
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

    tb = input$display_select_tbin
    sb = input$display_select_sbin

    flow_plot = summarise_flow_spatial(data, to_plot, tb, sb, metacell_types, metacell_names, from = input$from_to == 'From')
    plot_focus = flow_plot$plot_focus[,c('sbin2', 'flow')]
    plot_dynamics = flow_plot$plot_dynamics

    if(main){
        
        bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus, color)
        bp = flow_arrow_plot(bp, data, tb, plot_dynamics)
        bp = bp + gg_theme() + ggtitle(sprintf('%s \n time bin %s %s', to_plot, tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb))))

        return(bp)

    }else{

        p = list()
        idx = 1

        if(input$from_to == 'From'){

            to_flows = unique(plot_dynamics$ent2)
            to_flows = c(to_plot, to_flows[to_flows != to_plot])

            flow_frac = sum(plot_focus$flow)

            for(to_flow in to_flows){
                # "to flow" related to focus
                plot_focus_to_flow = plot_dynamics[plot_dynamics$ent2 == to_flow,] %>% group_by(sbin2) %>% summarise(flow = sum(flow))

                if(sum(plot_focus_to_flow$flow)/flow_frac >= f_th){
                    # browser()

                    # focus related to "to flow"
                    plot_focus = plot_dynamics[plot_dynamics$ent2 == to_flow,] %>% group_by(sbin1) %>% summarise(flow = sum(flow))

                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus, color)
                    bp = flow_arrow_plot(bp, data, tb, plot_dynamics[plot_dynamics$ent2 == to_flow,])
                    bp = bp + gg_theme() + 
                        ggtitle(sprintf('[%s]->%s \n time bin [%s] -> %s %s', 
                                to_plot, to_flow, tb, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb))))
                    p[[idx]] = bp

                    # "to flow" related to focus
                    flow_plot = summarise_flow_spatial(data, to_flow, as.character(as.numeric(tb)+1), sb, metacell_types, metacell_names, from = F)
                    plot_dynamics_to_flow = flow_plot$plot_dynamics
                    plot_dynamics_to_flow = plot_dynamics_to_flow[plot_dynamics_to_flow$ent1 == to_plot,]

                    if(input$mode == 'Types'){
                        cell_type_colors_tab = cell_type_colors()
                        color2 = cell_type_colors_tab[cell_type_colors_tab$cell_type == to_flow,]$color
                    }else if(input$mode == 'SMCs'){
                        metacell_types_tab = metacell_types()
                        color2 = metacell_types_tab[metacell_types_tab$metacell == to_flow,]$mc_col
                    }

                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus_to_flow, color2)
                    bp = flow_arrow_plot(bp, data, tb, plot_dynamics_to_flow)
                    bp = bp + gg_theme() +
                        ggtitle(sprintf('%s->[%s] \n time bin [%s] -> %s %s', 
                                    to_plot, to_flow, tb, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    p[[idx + 1]] = bp

                    # general to flow
                    plot_focus_to_flow_tot = flow_plot$plot_focus[,c('sbin2', 'flow')]
                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus_to_flow_tot, color2)
                    bp = bp + gg_theme() + 
                        ggtitle(sprintf('%s \n time bin %s %s', 
                                        to_flow, as.character(as.numeric(tb)+1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    p[[idx + 2]] = bp

                    idx  = idx + 3
                }
                    
            }
        }else{
            # browser()
            from_flows = unique(plot_dynamics$ent1)
            from_flows = c(to_plot, from_flows[from_flows != to_plot])

            flow_frac = sum(plot_focus$flow)

            for(from_flow in from_flows){
                # "to flow" related to focus
                plot_focus_from_flow = plot_dynamics[plot_dynamics$ent1 == from_flow,] %>% group_by(sbin1) %>% summarise(flow = sum(flow))

                if(sum(plot_focus_from_flow$flow)/flow_frac >= f_th){

                    # focus related to "to flow"
                    plot_focus = plot_dynamics[plot_dynamics$ent1 == from_flow,] %>% group_by(sbin2) %>% summarise(flow = sum(flow))

                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus, color)
                    bp = flow_arrow_plot(bp, data, tb, plot_dynamics[plot_dynamics$ent1 == from_flow,])
                    bp = bp + gg_theme() + 
                        ggtitle(sprintf('[%s]->%s \n time bin %s -> [%s] %s', 
                                from_flow, to_plot, as.character(as.numeric(tb)-1), tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb))))
                    p[[idx + 2]] = bp

                    # "to flow" related to focus
                    plot_focus_from_flow = plot_dynamics[plot_dynamics$ent1 == from_flow,] %>% group_by(sbin1) %>% summarise(flow = sum(flow))

                    if(input$mode == 'Types'){
                        cell_type_colors_tab = cell_type_colors()
                        color2 = cell_type_colors_tab[cell_type_colors_tab$cell_type == from_flow,]$color
                    }else if(input$mode == 'SMCs'){
                        metacell_types_tab = metacell_types()
                        color2 = metacell_types_tab[metacell_types_tab$metacell == from_flow,]$mc_col
                    }

                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus_from_flow, color2)
                    bp = flow_arrow_plot(bp, data, tb, plot_dynamics[plot_dynamics$ent1 == from_flow,])
                    bp = bp + gg_theme() +
                        ggtitle(sprintf('%s->[%s] \n time bin %s -> [%s] %s', 
                                    from_flow, to_plot, as.character(as.numeric(tb)-1), tb, ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    p[[idx + 1]] = bp

                    # general from flow
                    flow_plot = summarise_flow_spatial(data, from_flow, as.character(as.numeric(tb)-1), sb, metacell_types, metacell_names, from = F)
                    plot_focus_from_flow_tot = flow_plot$plot_focus[,c('sbin2', 'flow')]

                    bp = scatter_beatle_plot(input, data, tbin_time, tb, plot_focus_from_flow_tot, color2)
                    bp = bp + gg_theme() + 
                        ggtitle(sprintf('%s \n time bin %s %s', 
                                        from_flow, as.character(as.numeric(tb)-1), ifelse(sb == 'ALL', '', sprintf('\n %s', sb)))) 
                    p[[idx]] = bp

                    idx  = idx + 3
                }   
            }
        }
        # browser()
    
        g = arrangeGrob(grobs = lapply(p, ggplotGrob), ncol = 3,
                        top = sprintf('%s %s %s', to_plot, tb, ifelse(sb == 'ALL', '', sb)))

        return(grid.draw(g))
    }
    
}

scatter_beatle_plot = function(input, data, tbin_time, tb, to_plot_scatter, color, N = 10000){

    # browser()
    spat_model = round(data$tbin_time[[as.character(tb)]],1)
    proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))

    colnames(to_plot_scatter) = c('sbin', 'f')
    to_plot_scatter = to_plot_scatter %>% complete(sbin = unique(data$ed$sbin2), fill = list(f = 0))

    n_plot = to_plot_scatter$f
    names(n_plot) = to_plot_scatter$sbin
    n_plot[is.na(n_plot)] = 0
    n_plot = round(n_plot*N)

    poly_list <- split(proj$sbins, proj$sbins$sbin)
    
    points_list = data.frame(x = double(), y = double())
    points_df <- Map(generate_points_in_poly, poly_list, n_plot)
    points_df = do.call(rbind, points_df)

    p = proj$heatmap(fill = "white") + geom_point(data = points_df, aes(x = x, y = y), colour = color, alpha = 0.25)

    return(p)
}

runif_polygon <- function(n, poly, smooth=0, bounding=NULL, ncores=1)
{
    xmin <- min(poly$x)
    xmax <- max(poly$x)
    ymin <- min(poly$y)
    ymax <- max(poly$y)

    pts <- data.frame(x=runif(n*2, xmin, xmax), y=runif(n*2, ymin, ymax))
    pts <- pts[mingeom::within_pt_pol(poly, pts, ncores),]

    if (smooth > 0) {
        pts <- pts %>%
            mutate(x=x+rnorm(nrow(pts), 0, smooth)) %>%
            mutate(y=y+rnorm(nrow(pts), 0, smooth))

        if (!is.null(bounding)) {
            pts <- pts[mingeom::within_pt_pol(bounding, pts, ncores),]
        }
    }

    if (nrow(pts) < n) {
        pts <- bind_rows(pts, runif_polygon(n-nrow(pts), poly, smooth, bounding, ncores))
    }
    else {
        pts <- pts[1:n,]
    }

    return(pts)
}

generate_points_in_poly = function(df, n) {
    coords <- as.matrix(rbind(df[, c("x", "y")], df[1, c("x", "y")]))
    points = runif_polygon(n, as.data.frame(coords))
    return(points)            
}

flow_arrow_plot = function(p, data, tb, plot_dynamics, N = 10000){ 

        spat_model = round(data$tbin_time[[as.character(tb)]],1)
        proj <- spatula::eeSbinsModel$new(sprintf("s%s",sprintf("%.1f",spat_model)))
        
        sbin_centers = proj$sbins %>% group_by(sbin) %>% summarise(x = (max(x)+min(x))/2, y = (max(y)+min(y))/2)
        sbin_centers[sbin_centers$sbin == 'L1',]$x = sbin_centers[sbin_centers$sbin == 'L2',]$x - 0.02
        sbin_centers[sbin_centers$sbin == 'L6',]$x = sbin_centers[sbin_centers$sbin == 'L5',]$x + 0.02
        sbin_centers[sbin_centers$sbin == 'M6',]$x = sbin_centers[sbin_centers$sbin == 'L3',]$x - 0.02
        sbin_centers[sbin_centers$sbin == 'M7',]$x = sbin_centers[sbin_centers$sbin == 'L4',]$x + 0.02
        sbin_centers[sbin_centers$sbin == 'X1_1',]$y = 0.07
        sbin_centers[sbin_centers$sbin == 'X1_4',]$y = 0.07
        
        flow_dir = plot_dynamics %>% group_by(sbin1, sbin2) %>% summarise(flow = sum(flow))
        flow_dir = flow_dir[!flow_dir$sbin1 == flow_dir$sbin2 & flow_dir$flow > 1e-4,]

        x_from = c()
        x_to = c()
        y_from = c()
        y_to = c()

        if(nrow(flow_dir) > 0){
            for(row in seq(1, nrow(flow_dir))){

                x_from = c(x_from, sbin_centers[sbin_centers$sbin == flow_dir[row,]$sbin1, ]$x)
                y_from = c(y_from, sbin_centers[sbin_centers$sbin == flow_dir[row,]$sbin1, ]$y)
                x_to = c(x_to, sbin_centers[sbin_centers$sbin == flow_dir[row,]$sbin2, ]$x)
                y_to = c(y_to, sbin_centers[sbin_centers$sbin == flow_dir[row,]$sbin2, ]$y)

                dx = ifelse(x_to > x_from, 1, -1)
                dy = ifelse(y_to > y_from, 1, -1)

                # x_from = x_from + 0.001*dx
                # x_to = x_to - 0.005*dx
                # y_from = y_from + 0.001*dy
                # y_to = y_to - 0.005*dy

            }       
            # browser()
            color_scales = round(flow_dir$flow*N)
            colors_range = colorRampPalette(c("white", "blue"))(max(color_scales))

            p = p + geom_segment(aes(x =x_from, y = y_from, xend = x_to, yend = y_to, colour = color_scales), 
            arrow = arrow(length = unit(0.4, 'cm'), type = 'closed')) + theme(legend.position = "none")
            scale_colour_gradient(low = "white", high = "darkblue")
        }
        # ggsave(filename = 'test6.pdf', plot = p, device = "pdf")

        return(p)
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

    display_selector_tbin <- function(input, output, session, time_bins, ns) {
        output$display_select_tbin <- renderUI({
            tagList(
                shinyWidgets::pickerInput(ns("display_select_tbin"), "time bin",
                    selected = time_bins[2],
                    choices = time_bins,
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                )
            )
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