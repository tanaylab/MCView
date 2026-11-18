# Convert the legacy .qs cache files of MCView projects (imported with MCView < 0.2.54) to .qs2, which take about
# half the space and which MCView >= 0.2.56 reads faster (with several threads when qs2 has TBB).
#
#   Rscript dev/convert_qs_to_qs2.R [--dry-run] [--overwrite] [--remove-qs] <project or bundle dir> ...
#
# For each <name>.qs in <project>/cache/<dataset> and in its atlas/ subdirectory, reads it with qs, writes it with
# qs2 to a hidden temporary file in the same directory, reads that back (validating its checksum), checks that it
# is identical to the object read from the .qs (or serializes identically, see same()), and only then renames it to
# <name>.qs2. When both files exist, MCView reads the .qs2. The .qs is kept, the new .qs2 gets its permissions, and
# nothing else is written. Needs the qs and qs2 packages, not MCView. A run that is killed can leave a hidden
# .<name>.qs2-XXXX file, which MCView doesn't read: delete it.
#
#   --dry-run    list the .qs files that would be converted, with their sizes, without reading or writing anything
#   --overwrite  also convert the .qs files that have a .qs2 next to them, replacing it. Careful: a .qs2 that a
#                newer MCView wrote next to an old .qs (e.g. update_metadata()) is newer than the .qs.
#   --remove-qs  delete each .qs once a .qs2 identical to it is in place (converted now or found). Use only when no
#                app process of the bundle is running (e.g. right after restarting it): it may be reading the .qs.
#
# Skips bundles whose app.R loads its own MCView from a `code` directory: MCView < 0.2.54 can't read .qs2 and
# fails at startup when a .qs2 is next to the .qs files. `code` counts even when it is a symlink that doesn't
# resolve here (e.g. to /MCView/0.2.5, which exists in the server's container). Also skips dataset and atlas
# directories that are symlinks to outside the project, and directories it can't write to.
args <- commandArgs(TRUE)
opts <- c("--dry-run", "--overwrite", "--remove-qs")
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
remove_qs <- "--remove-qs" %in% args
paths <- setdiff(args, opts)
if (length(paths) == 0 || any(startsWith(paths, "--"))) {
    stop("usage: Rscript dev/convert_qs_to_qs2.R [--dry-run] [--overwrite] [--remove-qs] <project or bundle dir> ...")
}
if (!dry_run && !requireNamespace("qs", quietly = TRUE)) stop("the qs package is needed to read the .qs files")

threads <- 8
# qs2 built without TBB reads and writes with one thread and warns
qs2_call <- function(expr) {
    withCallingHandlers(expr, warning = function(w) if (grepl("TBB not available", conditionMessage(w))) invokeRestart("muffleWarning"))
}
fmt_size <- function(bytes) if (is.na(bytes) || bytes >= 1e9) sprintf("%.2f GB", bytes / 1e9) else sprintf("%.1f MB", bytes / 1e6)
# identical() compares environments by address, so it never matches two reads of an object that holds one (tables
# imported from anndata keep a dead reticulate object in a pandas.index attribute): compare their serializations then
same <- function(x, y) identical(x, y) || identical(serialize(x, NULL), serialize(y, NULL))

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(dirname(script), "own_mcview.R"))

# Convert one file. Returns what it did, or fails.
convert <- function(qs_file) {
    qs2_file <- paste0(qs_file, "2")
    keep_qs2 <- file.exists(qs2_file) && !overwrite
    if (keep_qs2 && !remove_qs) {
        qs2_call(qs2::qs_read(qs2_file, validate_checksum = TRUE, nthreads = threads))
        return("skipped: has a valid .qs2")
    }
    obj <- qs::qread(qs_file, nthreads = threads)
    if (keep_qs2) {
        if (!same(qs2_call(qs2::qs_read(qs2_file, validate_checksum = TRUE, nthreads = threads)), obj)) {
            stop("the .qs2 next to it differs from the .qs, kept both")
        }
        status <- "has an identical .qs2"
    } else {
        tmp <- tempfile(paste0(".", basename(qs2_file), "-"), tmpdir = dirname(qs_file))
        on.exit(unlink(tmp))
        qs2_call(qs2::qs_save(obj, tmp, nthreads = threads))
        if (!same(qs2_call(qs2::qs_read(tmp, validate_checksum = TRUE, nthreads = threads)), obj)) {
            stop("the .qs2 read back differs from the .qs")
        }
        Sys.chmod(tmp, file.mode(qs_file), use_umask = FALSE)
        if (!file.rename(tmp, qs2_file)) stop("could not rename the temporary file to ", qs2_file)
        status <- "converted"
    }
    if (remove_qs) {
        if (unlink(qs_file) != 0 || file.exists(qs_file)) stop("could not remove the .qs")
        status <- paste(status, "and removed the .qs")
    }
    status
}

convert_project <- function(path) {
    if (dir.exists(file.path(path, "project", "cache"))) {
        bundle <- path
        project <- file.path(path, "project")
    } else {
        bundle <- if (basename(path) == "project") dirname(path) else path
        project <- path
    }
    cat("\n== ", path, "\n", sep = "")
    why <- own_mcview(bundle)
    if (!is.null(why)) {
        cat("  skipped: ", why, "\n", sep = "")
        return(NULL)
    }
    cache <- file.path(project, "cache")
    if (!dir.exists(cache)) {
        cat("  FAILED: no cache directory\n")
        return(data.frame(bundle = path, failed = 1, todo = 0, todo_bytes = 0, unwritable = 0, unwritable_bytes = 0))
    }

    project_real <- normalizePath(project, winslash = "/")
    datasets <- list.files(cache, full.names = TRUE)
    dirs <- c(datasets, file.path(datasets, "atlas"))
    dirs <- sort(dirs[dir.exists(dirs)])
    res <- list()
    for (dir in dirs) {
        qs_files <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        if (length(qs_files) == 0) next
        rel <- substring(dir, nchar(cache) + 2)
        real <- normalizePath(dir, winslash = "/")
        if (!startsWith(real, paste0(project_real, "/"))) {
            cat("  ", rel, ": skipped, a symlink to ", real, ", outside the project (convert that project instead)\n", sep = "")
            next
        }
        if (file.access(dir, 3) != 0) { # can't create a file in it
            info <- file.info(dir)
            cat(sprintf(
                "  %s: skipped, can't create files in it (owner %s, mode %s): %d .qs files, %s\n",
                rel, info$uname, format(info$mode), length(qs_files), fmt_size(sum(file.size(qs_files)))
            ))
            res[[length(res) + 1]] <- data.frame(file = qs_files, qs = file.size(qs_files), qs2 = NA, seconds = 0, status = "skipped: directory not writable")
            next
        }
        for (qs_file in qs_files) {
            name <- file.path(rel, basename(qs_file))
            if (dry_run) {
                status <- if (file.exists(paste0(qs_file, "2")) && !overwrite) "has a .qs2" else "convert"
                cat(sprintf("  %-50s %10s  %s%s\n", name, fmt_size(file.size(qs_file)), status, if (remove_qs) ", then remove the .qs if identical" else ""))
                res[[length(res) + 1]] <- data.frame(file = qs_file, qs = file.size(qs_file), qs2 = NA, seconds = 0, status = status)
                next
            }
            size <- file.size(qs_file)
            seconds <- system.time(
                status <- tryCatch(convert(qs_file), error = function(e) paste("failed:", conditionMessage(e)))
            )[["elapsed"]]
            gc() # free the object before reading the next file
            qs2_size <- if (startsWith(status, "converted")) file.size(paste0(qs_file, "2")) else NA
            cat(sprintf("  %-50s %10s -> %10s %6.1fs  %s\n", name, fmt_size(size), if (is.na(qs2_size)) "" else fmt_size(qs2_size), seconds, status))
            res[[length(res) + 1]] <- data.frame(file = qs_file, qs = size, qs2 = qs2_size, seconds = seconds, status = status)
        }
    }
    if (length(res) == 0) {
        cat("  no .qs files\n")
        return(NULL)
    }
    res <- do.call(rbind, res)
    failed <- startsWith(res$status, "failed")
    if (dry_run) {
        todo <- res$status == "convert"
        cat(sprintf("  would convert %d of %d .qs files: %s\n", sum(todo), nrow(res), fmt_size(sum(res$qs[todo]))))
    } else {
        done <- startsWith(res$status, "converted")
        removed <- grepl("removed the .qs", res$status)
        cat(sprintf(
            "  converted %d files, %s -> %s, in %.0fs; skipped %d; removed %d .qs (%s); %d failed\n",
            sum(done), fmt_size(sum(res$qs[done])), fmt_size(sum(res$qs2[done])), sum(res$seconds), sum(startsWith(res$status, "skipped")),
            sum(removed), fmt_size(sum(res$qs[removed])), sum(failed)
        ))
    }
    unwritable <- res$status == "skipped: directory not writable"
    data.frame(
        bundle = path, failed = sum(failed), todo = sum(res$status == "convert"), todo_bytes = sum(res$qs[res$status == "convert"]),
        unwritable = sum(unwritable), unwritable_bytes = sum(res$qs[unwritable])
    )
}

summary <- do.call(rbind, lapply(paths, convert_project))
if (!is.null(summary)) {
    if (dry_run) {
        cat(sprintf("\nwould convert %d .qs files, %s, in %d projects\n", sum(summary$todo), fmt_size(sum(summary$todo_bytes)), sum(summary$todo > 0)))
    }
    if (any(summary$unwritable > 0)) {
        cat(sprintf("skipped %d .qs files, %s, in directories it can't write to\n", sum(summary$unwritable), fmt_size(sum(summary$unwritable_bytes))))
    }
    if (any(summary$failed > 0)) {
        cat("\nfailures in:", summary$bundle[summary$failed > 0], sep = "\n  ")
        cat("\n")
    }
}
quit(save = "no", status = if (is.null(summary) || all(summary$failed == 0)) 0 else 1)
