#!/usr/bin/env Rscript

#  option --all        draw all curves as jpg
#  option --only=name  draw only the curves of test "name" as jpg

args=commandArgs(trailingOnly=TRUE)
draw_all <- FALSE
only_one <- FALSE
if (length(args) > 0) {
    for (c_arg in 1:length(args)){
        if (args[c_arg] == "--all") {
            draw_all <- TRUE
        }
        if (startsWith(args[c_arg], '--only=')) {
            only_one <- TRUE
            long_name <- args[c_arg]
        }
    }
}

if (isTRUE(draw_all)) {
    filenames <- list.files("amplicon3", pattern="*.tmp$|*_output$", full.names=TRUE)
}
if (isTRUE(only_one)) {
    only_one <- FALSE
    ass_file_name <- paste("amplicon3/", substr(long_name, 8, 1000), sep = "",  collapse = "")
    if (file.exists(ass_file_name)) {
        only_one <- TRUE
        filenames <- list(ass_file_name)
    }
}

if (isTRUE(draw_all) | isTRUE(only_one)) {
    for (filename in filenames) {
        outname <- gsub(".tmp", "_tmp.jpg", filename)
        outname <- gsub("_output", "_output.jpg", outname)

        cat("Creating: ", outname, "\n")

        file_str <- readChar(filename, file.info(filename)$size)
        file_arr1 <- unlist(strsplit(file_str, "\n"))
        file_arr2 <- strsplit(file_arr1, "=")

        found_temp <- FALSE
        found_melt <- FALSE
        found_deri <- FALSE

        for (lin in 1:length(file_arr2)) {
            if (file_arr2[[lin]][[1]] == "AMPLICON_TEMPERATURES") {
                temp <- unlist(strsplit(file_arr2[[lin]][[2]], ","))
                found_temp <- TRUE
            }
            if (file_arr2[[lin]][[1]] == "AMPLICON_MELT_CURVE") {
                melt <- unlist(strsplit(file_arr2[[lin]][[2]], ","))
                found_melt <- TRUE
            }
            if (file_arr2[[lin]][[1]] == "AMPLICON_DERIVATIVE_CURVE") {
                deriv <- unlist(strsplit(file_arr2[[lin]][[2]], ","))
                found_deri <- TRUE
            }

        }

        if (found_temp & found_melt & found_deri) {
            jpeg(file=outname)
            plot(temp, melt, type = "l", col="grey")
            par(new = TRUE)
            plot(temp, deriv, type = "l", col="blue", axes = FALSE, bty = "n", xlab = "", ylab = "")
            axis(side=4, at = pretty(range(deriv)))
            mtext("z", side=4, line=3)
            dev.off()
        }
    }
}
