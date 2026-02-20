# Helpers for cleaning markdown content for display.

clean_about_markdown <- function(text) {
    if (is.null(text) || !nzchar(text)) {
        return(text)
    }

    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    if (length(lines) == 0) {
        return(text)
    }

    # Strip YAML front matter.
    if (trimws(lines[1]) == "---") {
        end_idx <- which(trimws(lines[-1]) == "---")
        if (length(end_idx) > 0) {
            end_idx <- end_idx[1] + 1
            if (end_idx < length(lines)) {
                lines <- lines[(end_idx + 1):length(lines)]
            } else {
                lines <- character()
            }
        }
    }

    # Strip knitr setup chunks (e.g. opts_chunk$set).
    remove_idx <- rep(FALSE, length(lines))
    i <- 1
    while (i <= length(lines)) {
        if (grepl("^```\\{r", lines[i])) {
            end_rel <- which(grepl("^```\\s*$", lines[(i + 1):length(lines)]))
            if (length(end_rel) == 0) {
                break
            }
            end_idx <- i + end_rel[1]
            chunk_lines <- lines[i:end_idx]
            if (any(grepl("knitr::opts_chunk\\$set", chunk_lines))) {
                remove_idx[i:end_idx] <- TRUE
            }
            i <- end_idx + 1
            next
        }
        i <- i + 1
    }

    lines <- lines[!remove_idx]
    while (length(lines) > 0 && !nzchar(trimws(lines[1]))) {
        lines <- lines[-1]
    }

    paste(lines, collapse = "\n")
}
