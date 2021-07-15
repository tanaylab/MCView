#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
args <- list()
if (length(argv) < 1) {
    cat("Usage: app.R <project> <port> <host>\n", sep = "", file = stderr())
    q(save = "no", status = 2)
}


args$project <- argv[1]

if (length(argv) > 1){
    args$port <- as.integer(argv[2])
}

if (length(argv) > 2){
    args$host <- argv[3]
}


MCView::run_app(project = args$project, port = args$port, host = args$host)