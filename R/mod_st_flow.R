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
            uiOutput(ns("norm_flow")),
            fileInput(ns("load_projection"),
                label = NULL,
                buttonLabel = "Load 2D layout",
                multiple = FALSE,
                accept =
                    c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        "text/tab-separated-values",
                        ".csv",
                        ".tsv"
                    )
            )
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

            output$Temporal_Flow = renderPlot({plot_temporal_flow_bars(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors)})
        }
    )
}

plot_temporal_flow_bars = function(input, output, session, dataset, data, metacell_types, metacell_names, cell_type_colors){
    
    req(input$mode)

    flow_to = data$flow_to
    flow_from = data$flow_from

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

        flow_from = flow_from %>% group_by(time_bin, ent1, ent2) %>% summarise(f = sum(f, na.rm = T))
        flow_from = flow_from %>% group_by(time_bin, ent2) %>% mutate(f_norm = f/sum(f))

        flow_to = flow_to %>% group_by(time_bin, ent1, ent2) %>% summarise(f = sum(f, na.rm = T))
        flow_to = flow_to %>% group_by(time_bin, ent1) %>% mutate(f_norm = f/sum(f))

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
   
    flow_to[is.na(flow_to$ent2),]$ent2 = 'sink'
    flow_to[is.na(flow_to$ent1),]$ent1 = 'source'

    subset_flow_from = flow_from[flow_from$ent1 == selected,]
    subset_flow_from = subset_flow_from[subset_flow_from$f > 0.01 | subset_flow_from$ent2 == selected,]

    subset_flow_to = flow_to[flow_to$ent2 == selected,]
    subset_flow_to = subset_flow_to[subset_flow_to$f > 0.01 | subset_flow_to$ent1 == selected,]

    subset_flow = data.frame(time_bin = c(subset_flow_from$time_bin, subset_flow_to$time_bin),
                             smc = c(subset_flow_from$ent2, subset_flow_to$ent1),
                             f = c(subset_flow_from$f, -subset_flow_to$f),
                             f_norm = c(subset_flow_from$f_norm, -subset_flow_to$f_norm))

    subset_flow$smc = factor(subset_flow$smc, levels = unique(subset_flow$smc))
    subset_flow$time_bin = factor(subset_flow$time_bin, levels = as.character(sort(as.numeric(unique(subset_flow$time_bin)))))
    subset_flow$selected = subset_flow$smc == selected

    if(input$norm_flow){
       subset_flow$flow_plot = subset_flow$f_norm
    }else if(!input$norm_flow){
       subset_flow$flow_plot = subset_flow$f
    }
    x_lim = max(abs(min(subset_flow$flow_plot, na.rm = T)), max(subset_flow$flow_plot, na.rm = T))
    g = ggplot(subset_flow, aes(y=smc, x=flow_plot, fill = smc)) + 
                geom_bar(stat = "identity", aes(color = selected), linewidth = 1.5) + 
                facet_wrap(~time_bin, ncol = 4) + 
                scale_fill_manual(values=ctype_color) + 
                scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"))+
                theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1),
                        strip.text = element_text()) +
                geom_vline(xintercept = 0, color = "black", linewidth = 1) +
                xlim(-x_lim, x_lim) + guides(color = "none")

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
         shinyWidgets::switchInput(
                inputId = ns("norm_flow"),
                label = "Normalize flow",
                size = "sm",
                value = TRUE
            )
        }
)}
