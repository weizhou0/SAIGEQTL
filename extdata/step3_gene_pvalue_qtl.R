#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors = F)
library(SAIGEQTL)
BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
print(sessionInfo())

option_list <- list(
  make_option("--assocFile",
    type = "character", default = "",
    help = "Path to file output by step 2 that contains p-values of single-variants."
  ),
  make_option("--weightFile",
    type = "character", default = "",
    help = "Path to file contains weights for each marker. The file has two columns: one with marker IDs that are used to match with MarkerID in the assocFile (output by step 2) and the other column contains weights used for each marker. The file needs to have a header with column names 'MarkerID' and 'weight'. If not specified, equal weights will be used for all markers. "
  ),
  make_option("--geneName",
    type = "character", default = "",
    help = "gene name that will be included in the output. If not specified, 'gene' will be used"
  ),
  make_option("--genePval_outputFile", type = "character", default = "", help = "Path to the output file containing gene p-value calcuated using ACAT test")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

if (BLASctl_installed) {
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}

if (!file.exists(opt$assocFile)) {
  stop("opt$assocFile ", opt$assocFile, "does not exist\n")
} else {
  data <- fread(opt$assocFile, header = T, select = c("MarkerID", "p.value"))
}
if (opt$weightFile != "") {
  if (!file.exists(opt$weightFile)) {
    stop("opt$weightFile ", opt$weightFile, "does not exist\n")
  } else {
    dataw <- fread(opt$weightFile, header = T)
    setkey(dataw, "MarkerID")
    setkey(data, "MarkerID")
    datam <- merge(data, dataw, by.x = "MarkerID", by.y = "MarkerID")
    datam <- datam[complete.cases(datam), ]
    cauchyPval <- get_CCT_pvalue(datam$p.value, datam$weight)
    top_pval <- min(datam$p.value)
    top_MarkerID <- datam$MarkerID[which.min(datam$p.value)]
  }
} else {
  print(head(data))
  data <- data[complete.cases(data), ]
  cauchyPval <- get_CCT_pvalue(data$p.value)
  top_pval <- min(data$p.value)
  top_MarkerID <- data$MarkerID[which.min(data$p.value)]
}

if (opt$geneName == "") {
  geneName <- "gene"
} else {
  geneName <- opt$geneName
}


cat("top_MarkerID ", top_MarkerID, "\n")
cat("top_pval ", top_pval, "\n")
resultdata <- data.frame(gene = geneName, ACAT_p = cauchyPval, top_MarkerID = top_MarkerID, top_pval = top_pval)



fwrite(resultdata, opt$genePval_outputFile, quote = F, sep = "\t", col.names = T, row.names = F)
