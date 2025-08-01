#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries

library(SAIGEQTL)

require(optparse) # install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile",
    type = "character", default = "",
    help = "Path to plink file for creating the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF and maxMissingRate can be used to specify the maximum missing rates  of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios"
  ),
  make_option("--bedFile",
    type = "character", default = "",
    help = "Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"
  ),
  make_option("--bimFile",
    type = "character", default = "",
    help = "Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"
  ),
  make_option("--famFile",
    type = "character", default = "",
    help = "Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"
  ),
  make_option("--phenoFile",
    type = "character", default = "",
    help = "Required. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarColList. All categorical covariates must also be included in covarColList."
  ),
  make_option("--phenoCol",
    type = "character", default = "",
    help = "Required. Column name for phenotype to be tested in the phenotype file, e.g CAD"
  ),
  make_option("--isRemoveZerosinPheno",
    type = "logical", default = FALSE,
    help = "Optional. Whether to remove zeros in the phenotype"
  ),
  make_option("--traitType", type = "character", default = "binary", help = "Required. binary or quantitative [default=binary]"),
  make_option("--invNormalize",
    type = "logical", default = FALSE,
    help = "Optional. Only for quantitative. Whether to perform the inverse normalization for the phenotype [default='FALSE']"
  ),
  make_option("--covarColList",
    type = "character", default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option("--sampleCovarColList",
    type = "character", default = "",
    help = "List of covariates that are on sample level (comma separated)"
  ),
  make_option("--longlCol",
    type = "character", default = "",
    help = ""
  ),
  make_option("--qCovarColList",
    type = "character", default = "",
    help = "List of categorical covariates (comma separated). All categorical covariates must also be in covarColList"
  ),
  make_option("--sampleIDColinphenoFile",
    type = "character", default = "IID",
    help = "Required. Column name of sample IDs in the phenotype file, e.g. IID"
  ),
  make_option("--cellIDColinphenoFile",
    type = "character", default = "",
    help = "Column name of cell IDs in the phenotype file, e.g. barcode"
  ),
  make_option("--tol",
    type = "numeric", default = 0.02,
    help = "Optional. Tolerance for fitting the null GLMM to converge [default=0.02]."
  ),
  make_option("--maxiter",
    type = "integer", default = 20,
    help = "Optional. Maximum number of iterations used to fit the null GLMM [default=20]."
  ),
  make_option("--tolPCG",
    type = "numeric", default = 1e-5,
    help = "Optional. Tolerance for PCG to converge [default=1e-5]."
  ),
  make_option("--maxiterPCG",
    type = "integer", default = 500,
    help = "Optional. Maximum number of iterations for PCG [default=500]."
  ),
  make_option("--nThreads",
    type = "integer", default = 1,
    help = "Optional. Number of threads (CPUs) to use [default=1]."
  ),
  make_option("--SPAcutoff",
    type = "numeric", default = 2,
    help = "Optional. Cutoff for the deviation of score test statistics from mean in the unit of sd to perform SPA [default=2]."
  ),
  make_option("--numRandomMarkerforVarianceRatio",
    type = "integer", default = 30,
    help = "Optional. An integer greater than 0. Number of markers to be randomly selected for estimating the variance ratio. The number will be automatically added by 10 until the coefficient of variantion (CV) for the variance ratio estimate is below ratioCVcutoff [default=30]."
  ),
  make_option("--skipModelFitting",
    type = "logical", default = FALSE,
    help = "Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"
  ),
  make_option("--skipVarianceRatioEstimation",
    type = "logical", default = FALSE,
    help = "Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"
  ),
  make_option("--memoryChunk",
    type = "numeric", default = 2,
    help = "Optional. Size (Gb) for each memory chunk [default=2]"
  ),
  make_option("--tauInit",
    type = "character", default = "0,0",
    help = "Optional. Initial values for tau. [default=0,0]"
  ),
  make_option("--LOCO",
    type = "logical", default = TRUE,
    help = "Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting the null model using the full GRM [default=TRUE]."
  ),
  make_option("--isLowMemLOCO",
    type = "logical", default = FALSE,
    help = "Whehter to output the model file by chromosome when LOCO=TRUE. If TRUE, the memory usage in Step 1 and Step 2 will be lower [default=FALSE]"
  ),
  make_option("--traceCVcutoff",
    type = "numeric", default = 0.0025,
    help = "Optional. Threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold [default=0.0025]."
  ),
  make_option("--nrun",
    type = "numeric", default = 30,
    help = "Number of rums in trace estimation. [default=30]"
  ),
  make_option("--ratioCVcutoff",
    type = "numeric", default = 0.001,
    help = "Optional. Threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold [default=0.001]"
  ),
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "Required. Path and prefix of the output files [default='~/']"
  ),
  make_option("--outputPrefix_varRatio",
    type = "character", default = "",
    help = "Optional. Path and prefix of the output the variance ratio file. if not specified, it will be the same as the outputPrefix"
  ),
  make_option("--IsOverwriteVarianceRatioFile",
    type = "logical", default = FALSE,
    help = "Optional. Whether to overwrite the variance ratio file if the file exist.[default='FALSE']"
  ),
  make_option("--sparseGRMFile",
    type = "character", default = "",
    help = "Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"
  ),
  make_option("--sparseGRMSampleIDFile",
    type = "character", default = "",
    help = "Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]"
  ),
  make_option("--isCateVarianceRatio",
    type = "logical", default = FALSE,
    help = "Required. Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE [default=FALSE]"
  ),
  make_option("--relatednessCutoff",
    type = "numeric", default = 0,
    help = "Optional. Threshold (minimum relatedness coefficient) to treat two samples as unrelated when the sparse GRM is used [default=0]"
  ),
  make_option("--cateVarRatioMinMACVecExclude",
    type = "character", default = "10,20.5",
    help = "Optional. vector of float. Lower bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='10,20.5']"
  ),
  make_option("--cateVarRatioMaxMACVecInclude",
    type = "character", default = "20.5",
    help = "Optional. vector of float. Higher bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='20.5']"
  ),
  make_option("--isCovariateTransform",
    type = "logical", default = TRUE,
    help = "Optional. Whether use qr transformation on covariates [default='TRUE']."
  ),
  make_option("--isDiagofKinSetAsOne",
    type = "logical", default = FALSE,
    help = "Optional. Whether to set the diagnal elements in GRM to be 1 [default='FALSE']."
  ),
  make_option("--useSparseGRMtoFitNULL",
    type = "logical", default = FALSE,
    help = "Optional. Whether to use sparse GRM to fit the null model [default='FALSE']."
  ),
  make_option("--useSparseGRMforVarRatio",
    type = "logical", default = FALSE,
    help = "Optional. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE"
  ),
  make_option("--minMAFforGRM",
    type = "numeric", default = 0.01,
    help = "Optional. Minimum MAF of markers used for GRM"
  ),
  make_option("--maxMissingRateforGRM",
    type = "numeric", default = 0.15,
    help = "Optional. Maximum missing rate of markers used for GRM"
  ),
  make_option("--minCovariateCount",
    type = "numeric", default = -1,
    help = "Optional. Binary covariates with a count less than minCovariateCount will be excluded from the model to avoid convergence issues [default=-1] (no covariates will be excluded)."
  ),
  make_option("--includeNonautoMarkersforVarRatio",
    type = "logical", default = FALSE,
    help = "Optional. Whether to allow for non-autosomal markers for variance ratio. [default, 'FALSE']"
  ),
  make_option("--FemaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified"
  ),
  make_option("--MaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified"
  ),
  make_option("--sexCol",
    type = "character", default = "",
    help = "Optional. Column name for sex in the phenotype file, e.g Sex"
  ),
  make_option("--FemaleCode",
    type = "character", default = "1",
    help = "Optional. Values in the column for sex in the phenotype file are used for females [default, '1']"
  ),
  make_option("--MaleCode",
    type = "character", default = "0",
    help = "Optional. Values in the column for sex in the phenotype file are used for males [default, '0']"
  ),
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"
  ),
  make_option("--SampleIDIncludeFile",
    type = "character", default = "",
    help = "Path to the file that contains one column for IDs of samples who will be include for null model fitting."
  ),
  make_option("--VmatFilelist",
    type = "character", default = "",
    help = "List of additional V (comma separated)"
  ),
  make_option("--VmatSampleFilelist",
    type = "character", default = "",
    help = "List of additional V (comma separated)"
  ),
  make_option("--useGRMtoFitNULL", type = "logical", default = TRUE, help = ""),
  make_option("--offsetCol",
    type = "character", default = NULL,
    help = "offset column"
  ),
  make_option("--varWeightsCol",
    type = "character", default = NULL,
    help = "variance weight column"
  ),
  make_option("--isStoreSigma",
    type = "logical", default = TRUE,
    help = "Optional. Whether to store the inv Sigma matrix. [default, 'TRUE']"
  ),
  make_option("--isShrinkModelOutput",
    type = "logical", default = TRUE,
    help = "Optional. Whether to remove unnecessary objects for step2 from the model output. [default, 'TRUE']"
  ),
  make_option("--isExportResiduals",
    type = "logical", default = FALSE,
    help = "Optional. Whether to export residual vector. [default, 'FALSE']"
  )
)


## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList, ",")[[1]]
qcovars <- strsplit(opt$qCovarColList, ",")[[1]]
scovars <- strsplit(opt$sampleCovarColList, ",")[[1]]
convertoNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) {
    return(NULL)
  })
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  } else {
    cat(stringOutput, " is ", y, "\n")
  }
  return(y)
}

tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
cateVarRatioMinMACVecExclude <- convertoNumeric(x = strsplit(opt$cateVarRatioMinMACVecExclude, ",")[[1]], "cateVarRatioMinMACVecExclude")
cateVarRatioMaxMACVecInclude <- convertoNumeric(x = strsplit(opt$cateVarRatioMaxMACVecInclude, ",")[[1]], "cateVarRatioMaxMACVecInclude")

BLASctl_installed <- require(RhpcBLASctl)
if (BLASctl_installed) {
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}


# set seed
set.seed(1)
fitNULLGLMM_multiV(
  plinkFile = opt$plinkFile,
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
  sparseGRMFile = opt$sparseGRMFile,
  sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
  phenoFile = opt$phenoFile,
  phenoCol = opt$phenoCol,
  isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
  sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
  cellIDColinphenoFile = opt$cellIDColinphenoFile,
  traitType = opt$traitType,
  outputPrefix = opt$outputPrefix,
  isCovariateOffset = opt$isCovariateOffset,
  nThreads = opt$nThreads,
  useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
  invNormalize = opt$invNormalize,
  covarColList = covars,
  qCovarCol = qcovars,
  tol = opt$tol,
  maxiter = opt$maxiter,
  tolPCG = opt$tolPCG,
  maxiterPCG = opt$maxiterPCG,
  SPAcutoff = opt$SPAcutoff,
  numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
  skipModelFitting = opt$skipModelFitting,
  skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
  memoryChunk = opt$memoryChunk,
  tauInit = tauInit,
  LOCO = opt$LOCO,
  isLowMemLOCO = opt$isLowMemLOCO,
  traceCVcutoff = opt$traceCVcutoff,
  nrun = opt$nrun,
  ratioCVcutoff = opt$ratioCVcutoff,
  outputPrefix_varRatio = opt$outputPrefix_varRatio,
  IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
  relatednessCutoff = opt$relatednessCutoff,
  isCateVarianceRatio = opt$isCateVarianceRatio,
  cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
  cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
  isCovariateTransform = opt$isCovariateTransform,
  isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
  minMAFforGRM = opt$minMAFforGRM,
  maxMissingRateforGRM = opt$maxMissingRateforGRM,
  minCovariateCount = opt$minCovariateCount,
  includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
  sexCol = opt$sexCol,
  FemaleCode = opt$FemaleCode,
  FemaleOnly = opt$FemaleOnly,
  MaleCode = opt$MaleCode,
  MaleOnly = opt$MaleOnly,
  SampleIDIncludeFile = opt$SampleIDIncludeFile,
  VmatFilelist = opt$VmatFilelist,
  VmatSampleFilelist = opt$VmatSampleFilelist,
  longlCol = opt$longlCol,
  useGRMtoFitNULL = opt$useGRMtoFitNULL,
  offsetCol = opt$offsetCol,
  varWeightsCol = opt$varWeightsCol,
  sampleCovarCol = scovars,
  isStoreSigma = opt$isStoreSigma,
  isShrinkModelOutput = opt$isShrinkModelOutput,
  isExportResiduals = opt$isExportResiduals
)

if (!opt$isCovariateOffset) {
  my_env <- new.env()
  load(paste0(opt$outputPrefix, ".rda"), envir = my_env)
  modglmm <- my_env$modglmm
  print(modglmm$theta)
  if (sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 1) {
    cat("All variance component parameter estiamtes are out of bounds, now try including all covariates as offset\n")
    opt$isCovariateOffset <- TRUE
    set.seed(1)
    fitNULLGLMM_multiV(
      plinkFile = opt$plinkFile,
      bedFile = opt$bedFile,
      bimFile = opt$bimFile,
      famFile = opt$famFile,
      useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
      sparseGRMFile = opt$sparseGRMFile,
      sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
      phenoFile = opt$phenoFile,
      phenoCol = opt$phenoCol,
      isRemoveZerosinPheno = opt$isRemoveZerosinPheno,
      sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
      cellIDColinphenoFile = opt$cellIDColinphenoFile,
      traitType = opt$traitType,
      outputPrefix = paste0(opt$outputPrefix, ".offset"),
      isCovariateOffset = opt$isCovariateOffset,
      nThreads = opt$nThreads,
      useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
      invNormalize = opt$invNormalize,
      covarColList = covars,
      qCovarCol = qcovars,
      tol = opt$tol,
      maxiter = opt$maxiter,
      tolPCG = opt$tolPCG,
      maxiterPCG = opt$maxiterPCG,
      SPAcutoff = opt$SPAcutoff,
      numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
      skipModelFitting = opt$skipModelFitting,
      skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
      memoryChunk = opt$memoryChunk,
      tauInit = tauInit,
      LOCO = opt$LOCO,
      isLowMemLOCO = opt$isLowMemLOCO,
      traceCVcutoff = opt$traceCVcutoff,
      nrun = opt$nrun,
      ratioCVcutoff = opt$ratioCVcutoff,
      outputPrefix_varRatio = opt$outputPrefix_varRatio,
      IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
      relatednessCutoff = opt$relatednessCutoff,
      isCateVarianceRatio = opt$isCateVarianceRatio,
      cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      isCovariateTransform = opt$isCovariateTransform,
      isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
      minMAFforGRM = opt$minMAFforGRM,
      maxMissingRateforGRM = opt$maxMissingRateforGRM,
      minCovariateCount = opt$minCovariateCount,
      includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
      sexCol = opt$sexCol,
      FemaleCode = opt$FemaleCode,
      FemaleOnly = opt$FemaleOnly,
      MaleCode = opt$MaleCode,
      MaleOnly = opt$MaleOnly,
      SampleIDIncludeFile = opt$SampleIDIncludeFile,
      VmatFilelist = opt$VmatFilelist,
      VmatSampleFilelist = opt$VmatSampleFilelist,
      longlCol = opt$longlCol,
      useGRMtoFitNULL = opt$useGRMtoFitNULL,
      offsetCol = opt$offsetCol,
      varWeightsCol = opt$varWeightsCol,
      sampleCovarCol = scovars,
      isStoreSigma = opt$isStoreSigma,
      isShrinkModelOutput = opt$isShrinkModelOutput,
      isExportResiduals = opt$isExportResiduals
    )
    my_env <- new.env()
    load(paste0(opt$outputPrefix, ".offset.rda"), envir = my_env)
    modglmm <- my_env$modglmm
    print(modglmm$theta)
    if (sum(modglmm$theta[2:length(modglmm$theta)]) <= 0 || sum(modglmm$theta[2:length(modglmm$theta)]) > 1) {
      cat("All variance component parameter estiamtes are out of bounds.\n")
      file.remove(paste0(opt$outputPrefix, ".offset.rda"))
      if (file.exists(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))) {
        file.remove(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))
        # Delete file if it exists
      } else {
        if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
          file.remove(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))
        }
      }
    } else {
      file.rename(paste0(opt$outputPrefix, ".offset.rda"), paste0(opt$outputPrefix, ".rda"))
      if (file.exists(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"))) {
        file.rename(paste0(opt$outputPrefix, ".offset.varianceRatio.txt"), paste0(opt$outputPrefix, ".varianceRatio.txt"))
        # Delete file if it exists
      } else {
        if (file.exists(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"))) {
          file.rename(paste0(opt$outputPrefix_varRatio, ".offset.varianceRatio.txt"), paste0(opt$outputPrefix_varRatio, ".varianceRatio.txt"))
        }
      }
    }
  }
}


if (BLASctl_installed) {
  # Restore originally configured BLAS thread count
  blas_set_num_threads(original_num_threads)
}
