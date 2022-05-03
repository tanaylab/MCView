get_md_attribute <- function(dataset, md, attr, default, atlas = FALSE) {
    metadata_colors <- get_mc_data(dataset, "metadata_colors", atlas = atlas)
    if (has_name(metadata_colors, md)) {
        if (is.null(attr)) {
            md_attr <- metadata_colors[[md]]
        } else {
            md_attr <- metadata_colors[[md]][[attr]]
        }

        if (!is.null(md_attr)) {
            return(md_attr)
        }
    }

    return(default)
}


get_metadata_colors <- function(dataset, md, colors = NULL, color_breaks = NULL, metadata = NULL, default_colors = c("white", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"), atlas = FALSE) {
    metadata <- metadata %||% get_mc_data(dataset, "metadata", atlas = atlas)
    if (md == "Clipboard") {
        return(c("selected" = "darkred", "unslected" = "lightgrey"))
    }
    if (is_numeric_field(metadata, md)) {
        colors <- colors %||% get_md_attribute(dataset, md, "colors", default_colors, atlas = atlas)
        color_breaks <- color_breaks %||% get_md_attribute(dataset, md, "breaks", NULL, atlas = atlas)

        if (is.null(color_breaks)) {
            min_val <- min(metadata[[md]], na.rm = TRUE)
            max_val <- max(metadata[[md]], na.rm = TRUE)

            if (min_val == max_val) {
                min_val <- min_val - 1e-5
            }

            color_breaks <- seq(min_val, max_val, length.out = length(colors))
        }

        return(
            list(colors = colors, breaks = color_breaks)
        )
    } else {
        colors <- colors %||% get_md_attribute(dataset, md, NULL, NULL, atlas = atlas)
        if (is.null(colors)) {
            categories <- unique(metadata[[md]])
            colors <- chameleon::distinct_colors(length(categories))$name
            names(colors) <- categories
        }
        return(colors)
    }
}
