#!/usr/bin/env Rscript

#  option -p  draw curves as jpg

args=commandArgs(trailingOnly=TRUE)
draw_pics <- FALSE
if (length(args) > 0) {
    for (c_arg in 1:length(args)){
        cat(c_arg, "\n")
        if (args[c_arg] == "-p") {
            draw_pics <- TRUE
        }
    }
}

cat("Installing Decipher library...\n")

if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
BiocManager::install("DECIPHER")

cat("Loading Decipher library...\n")
library("DECIPHER")

cat("Loading Amplicons...\n")
allSeqs <- read.table(file = 'amplicon3/amplicons.csv', header = FALSE, sep=";", stringsAsFactors=FALSE)

sName <- character(nrow(allSeqs) - 1)
rTemp <- numeric(nrow(allSeqs) - 1)

for (wCol in 1:nrow(allSeqs)){
    cat(paste("Processing ", allSeqs[wCol, 1], "...\n", sep = "",  collapse = ""))
    seq = allSeqs[wCol, 2]
    amp = DNAStringSet(seq)
    temp = (600:950) / 10

    prob = MeltDNA(amp,type = "positional probabilities",temps = temp, ions = 0.150)
    probTab = prob[[1]]
    head_nr <- rep(1:ncol(probTab))
    head_tab <- rbind(c(head_nr), probTab)
    col_temp <- append(temp, c(0), 0)
    final_tab <- cbind(c(col_temp), head_tab)
    file_name <- paste("amplicon3/decipher_", allSeqs[wCol, 1], "_heli.csv", sep = "",  collapse = "")
    write.table(final_tab, file = file_name, sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)

    melt = MeltDNA(amp,type = "melt curves",temps = temp, ions = 0.150)
    melt_tab <- cbind(c(temp), melt)
    melt_name <- paste("amplicon3/decipher_", allSeqs[wCol, 1], "_melt.csv", sep = "",  collapse = "")
    write.table(melt_tab, file = melt_name, sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)

    deriv = MeltDNA(amp,type = "derivative curves",temps = temp, ions = 0.150)
    deriv_tab <- cbind(c(temp), deriv)
    deriv_name <- paste("amplicon3/decipher_", allSeqs[wCol, 1], "_deriv.csv", sep = "",  collapse = "")
    write.table(deriv_tab, file = deriv_name, sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)

    if (isTRUE(draw_pics)) {
        jpeg(file=paste("amplicon3/decipher_", allSeqs[wCol, 1], "_curve.jpg", sep = "",  collapse = ""))
        plot(temp, melt, type = "l", col="grey")
        par(new = TRUE)
        plot(temp, deriv, type = "l", col="blue", axes = FALSE, bty = "n", xlab = "", ylab = "")
        axis(side=4, at = pretty(range(deriv)))
        mtext("z", side=4, line=3)
        dev.off()
    }

    sName[wCol - 1] <- allSeqs[wCol, 1]
    calcTemp <- temp[which.max(deriv)]
    rTemp[wCol - 1] <- calcTemp
    cat(paste(allSeqs[wCol, 1], ": ", calcTemp, "°C\n", sep = "",  collapse = ""))
}

finalData <- data.frame(sName, rTemp, stringsAsFactors=FALSE)
write.csv(finalData, "amplicon3/decipher_results.csv", row.names = FALSE)
