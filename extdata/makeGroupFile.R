#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries
# library(SAIGE, lib.loc="../../install_0.93")
# library(SAIGE, lib.loc="/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/installs_test")
library(SAIGEQTL)
require(optparse) # install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments

option_list <- list(
  make_option("--vcfFile",
    type = "character", default = "",
    help = "Path to vcf file."
  ),
  make_option("--vcfFileIndex",
    type = "character", default = "",
    help = "Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .csi file by tabix -p vcf csi file.vcf.gz"
  ),
  make_option("--vcfField",
    type = "character", default = "DS",
    help = "DS or GT, [default=DS]"
  ),
  make_option("--savFile",
    type = "character", default = "",
    help = "Path to the sav file."
  ),
  make_option("--savFileIndex",
    type = "character", default = "",
    help = "Path to the .s1r file (index of the sav file)."
  ),
  make_option("--bgenFile",
    type = "character", default = "",
    help = "Path to bgen file. Path to bgen file. Currently version 1.2 with 8 bit compression is supported"
  ),
  make_option("--bgenFileIndex",
    type = "character", default = "",
    help = "Path to the .bgi file (index of the bgen file)"
  ),
  make_option("--sampleFile",
    type = "character", default = "",
    help = "Path to the file that contains one column for IDs of samples in the dosage file. For version >= 0.38, this file is only needed for bgen files. "
  ),
  make_option("--bedFile",
    type = "character", default = "",
    help = "Path to bed file (PLINK)"
  ),
  make_option("--bimFile",
    type = "character", default = "",
    help = "Path to bim file (PLINK)"
  ),
  make_option("--famFile",
    type = "character", default = "",
    help = "Path to fam file (PLINK)"
  ),
  make_option("--AlleleOrder",
    type = "character", default = "alt-first",
    help = "alt-first or ref-first for bgen or PLINK files"
  ),
  make_option("--regionFile",
    type = "character", default = "",
    help = "Path to a file containing genome regions to extract genetic markers for making group files"
  ),
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "path and prefix to the output files [default='~/']"
  )
)


## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

set.seed(1)

# if(opt$plinkFile != ""){
#        bimFile = paste0(opt$plinkFile, ".bim")
#        bedFile = paste0(opt$plinkFile, ".bed")
#        famFile = paste0(opt$plinkFile, ".fam")
# }else{
# 	bimFile = opt$bimFile
# 	bedFile = opt$bedFile
# 	famFile = opt$famFile
# }



makeGroupFileforRegions(
  vcfFile = opt$vcfFile,
  vcfFileIndex = opt$vcfFileIndex,
  vcfField = opt$vcfField,
  savFile = opt$savFile,
  savFileIndex = opt$savFileIndex,
  bgenFile = opt$bgenFile,
  bgenFileIndex = opt$bgenFileIndex,
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  AlleleOrder = opt$AlleleOrder,
  regionFile = opt$regionFile,
  outputPrefix = opt$outputPrefix
)
