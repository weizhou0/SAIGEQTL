#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarCol. All categorical covariates must also be included in covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex").
#' @param eCovarCol vector of characters. Covariates of environmental factors/cell context to be used in the null model. All covariates listed in eCovarCol must be also in covarColList,  e,g c("cellType").
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file e.g. "IID".
#' @param cellIDColinphenoFile character. Column name for the cell IDs in the phenotype file e.g. "barcode".
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500.
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkersForVarRatio integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Initial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, TRUE
#' @param traceCVcutoff numeric. The threshold for coefficient of variation (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variation (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkersForVarRatio will be increased by 10. By default, 0.001
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix. By default, ""
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, ""
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, ""
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM.
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param useSparseGRMforVarRatio logical. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1'
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, ''
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, FALSE.
#' @param isStoreSigma logical. Whether to store sigma matrix. By default, FALSE. If number of individuals is greater than 10,000, this option may use large memory
#' @param isShrinkModelOutput logical. remove unnecessary objects for step2 from the model output. By default, FALSE.
#' @param isExportResiduals logical. export a residual vector. By default, FALSE.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV <- function(plinkFile = "",
                               bedFile = "",
                               bimFile = "",
                               famFile = "",
                               phenoFile = "",
                               phenoCol = "",
                               isRemoveZerosinPheno = FALSE,
                               traitType = "binary",
                               invNormalize = FALSE,
                               covarColList = NULL,
                               qCovarCol = NULL,
                               eCovarCol = NULL,
                               sampleCovarCol = NULL,
                               offsetCol = NULL,
                               varWeightsCol = NULL,
                               longlCol = "",
                               sampleIDColinphenoFile = "",
                               cellIDColinphenoFile = "",
                               tol = 0.02,
                               maxiter = 20,
                               tolPCG = 1e-5,
                               maxiterPCG = 500,
                               nThreads = 1,
                               SPAcutoff = 2,
                               numMarkersForVarRatio = 30,
                               skipModelFitting = FALSE,
                               memoryChunk = 2,
                               tauInit = c(0, 0),
                               LOCO = TRUE,
                               isLowMemLOCO = FALSE,
                               traceCVcutoff = 0.0025,
                               ratioCVcutoff = 0.001,
                               outputPrefix = "",
                               outputPrefix_varRatio = "",
                               IsOverwriteVarianceRatioFile = FALSE,
                               sparseGRMFile = "",
                               sparseGRMSampleIDFile = "",
                               numRandomMarkerforSparseKin = 1000,
                               relatednessCutoff = 0.125,
                               isCateVarianceRatio = FALSE,
                               cateVarRatioIndexVec = NULL,
                               cateVarRatioMinMACVecExclude = c(10, 20.5),
                               cateVarRatioMaxMACVecInclude = c(20.5),
                               isCovariateTransform = TRUE,
                               isDiagofKinSetAsOne = FALSE,
                               minCovariateCount = -1,
                               minMAFforGRM = 0.01,
                               maxMissingRateforGRM = 0.15,
                               useSparseGRMtoFitNULL = FALSE,
                               useSparseGRMforVarRatio = FALSE,
                               includeNonautoMarkersforVarRatio = FALSE,
                               sexCol = "",
                               FemaleCode = 1,
                               FemaleOnly = FALSE,
                               MaleCode = 0,
                               MaleOnly = FALSE,
                               SampleIDIncludeFile = "",
                               isCovariateOffset = FALSE,
                               skipVarianceRatioEstimation = FALSE,
                               nrun = 30,
                               VmatFilelist = "",
                               VmatSampleFilelist = "",
                               VcellmatFilelist = "",
                               VcellmatSampleFilelist = "",
                               useGRMtoFitNULL = TRUE,
                               isStoreSigma = FALSE,
                               isShrinkModelOutput = FALSE,
                               isExportResiduals = FALSE) {
  ## set up output files
  modelOut <- paste0(outputPrefix, ".rda")

  if (skipModelFitting) {
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
  } else {
    if (LOCO & isLowMemLOCO) {
      modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
    }
    file.create(modelOut, showWarnings = TRUE)
  }

  if (plinkFile != "") {
    bimFile <- paste0(plinkFile, ".bim")
    bedFile <- paste0(plinkFile, ".bed")
    famFile <- paste0(plinkFile, ".fam")
  }
  setgenoNULL()

  if (!useGRMtoFitNULL) {
    useSparseGRMtoFitNULL <- FALSE
    useSparseGRMforVarRatio <- FALSE
    LOCO <- FALSE
    # nThreads <- 1
    cat("No GRM will be used to fit the NULL model and nThreads is set to 1\n")
  }


  if (useSparseGRMtoFitNULL & bedFile == "") {
    cat("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated\n")
    skipVarianceRatioEstimation <- TRUE
  }

  if (!skipVarianceRatioEstimation) {
    SPAGMMATOut <- paste0(outputPrefix, "_", numMarkersForVarRatio, "markers.SAIGE.results.txt")
    # Check_OutputFile_Create(SPAGMMATOut)

    if (outputPrefix_varRatio == "") {
      outputPrefix_varRatio <- outputPrefix
    }
    varRatioFile <- paste0(outputPrefix_varRatio, ".varianceRatio.txt")

    if (!file.exists(varRatioFile)) {
      file.create(varRatioFile, showWarnings = TRUE)
    } else {
      if (!IsOverwriteVarianceRatioFile) {
        stop(
          "WARNING: The variance ratio file ", varRatioFile,
          " already exists. The new variance ratios will be output to ",
          varRatioFile, ". In order to avoid overwriting the file, please remove the ",
          varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify --IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)\n"
        )
      } else {
        cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
      }
    }
  } else {
    cat("Variance ratio estimation will be skipped\n.")
    useSparseGRMforVarRatio <- FALSE
  }


  if (useSparseGRMtoFitNULL) {
    # useSparseGRMforVarRatio = FALSE
    LOCO <- FALSE
    # nThreads <- 1
    # if (bedFile != "") {
    #  cat("sparse GRM will be used to fit the NULL model and nThreads is set to 1\n")
    # }
    cat("Leave-one-chromosome-out is not applied\n")
  }



  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!file.exists(sparseGRMFile)) {
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
    if (!file.exists(sparseGRMSampleIDFile)) {
      stop(
        "sparseGRMSampleIDFile ", sparseGRMSampleIDFile,
        " does not exist!"
      )
    }
  }


  if (nThreads > 1) {
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads will be used ", "\n")
  }
  set_g_omp_num_threads(nThreads)

  if (FemaleOnly & MaleOnly) {
    stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
  }

  if (FemaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_FemaleOnly")
    cat(
      "Female-specific model will be fitted. Samples coded as ",
      FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  } else if (MaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_MaleOnly")
    cat(
      "Male-specific model will be fitted. Samples coded as ",
      MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  }


  sampleListwithGeno <- NULL
  if ((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | !skipVarianceRatioEstimation) {
    if (!file.exists(bedFile)) {
      stop("ERROR! bed file does not exsit\n")
    }
    if (!file.exists(bimFile)) {
      stop("ERROR! bim file does not exsit\n")
    } else {
      if (LOCO) {
        chrVec <- data.table::fread(bimFile, header = F, data.table = F, select = 1)
        updatechrList <- updateChrStartEndIndexVec(chrVec)
        LOCO <- updatechrList$LOCO
        chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
        chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
      }
      if (!LOCO) {
        chromosomeStartIndexVec <- rep(NA, 22)
        chromosomeEndIndexVec <- rep(NA, 22)
      }
    }


    if (!file.exists(famFile)) {
      stop("ERROR! fam file does not exsit\n")
    } else {
      sampleListwithGenov0 <- data.table::fread(famFile, header = F, colClasses = list(character = 1:4), data.table = FALSE)
      colnames(sampleListwithGenov0) <- c(
        "FIDgeno", "IIDgeno",
        "father", "mother", "sex", "phe"
      )
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples have genotypes\n")
    }
  } else {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sampleListwithGenov0 <- data.table::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      sampleListwithGeno <- NULL
      sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
      sampleListwithGeno <- data.frame(sampleListwithGeno)
      sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
        by = 1
      )
      cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
    }
  }
  if (!file.exists(phenoFile)) {
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  } else {
    if (longlCol == "") {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile)
    } else {
      checkColList <- c(phenoCol, covarColList, sampleIDColinphenoFile, longlCol)
    }

    if (cellIDColinphenoFile != "") {
      cat(cellIDColinphenoFile, "is the cell ID column\n")
      checkColList <- c(checkColList, cellIDColinphenoFile)
    }

    if (length(offsetCol) > 0) {
      cat(offsetCol, "is the offset term\n")
      checkColList <- c(checkColList, offsetCol)
    }

    if (length(varWeightsCol) > 0) {
      cat(varWeightsCol, " is the weights for variance\n")
      checkColList <- c(checkColList, varWeightsCol)
    }


    ## check whether the phenotype file is large
    cmd <- paste0("du ", phenoFile, "| awk '{print $1}' > ", outputPrefix, "_", phenoCol, "_size_temp")
    system(cmd)
    datasize <- data.table::fread(paste0(outputPrefix, "_", phenoCol, "_size_temp"), header = F, data.table = F)
    isphenoFileLarge <- FALSE
    if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
      if (datasize[1, 1] > 200000) {
        isphenoFileLarge <- TRUE
      }
    } else {
      if (datasize[1, 1] > 500000) {
        isphenoFileLarge <- TRUE
      }
    }

    if (isphenoFileLarge) {
      catcmd <- ifelse(grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile), "gunzip -c ", "cat ")
      cmd <- paste0(catcmd, phenoFile, " | head -n 1 | sed 's/[\\t ]/\\n/g' | awk '{print $1\"\\t\"NR}' > ", outputPrefix, "_", phenoCol, "_lineNum_temp")
      system(cmd)

      checkColListDataFrame <- data.frame(colna = checkColList)
      phenoFilephenoCol_lineNum <- data.table::fread(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"), header = F, data.table = F)

      phenoFilephenoCol_lineNum_checkColList <- merge(checkColListDataFrame, phenoFilephenoCol_lineNum, by.x = 1, by.y = 1)

      write.table(phenoFilephenoCol_lineNum_checkColList[, 2], paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"), quote = F, col.names = F, row.names = F)

      cmdb <- paste0(catcmd, phenoFile, " | cut -f $(tr '\\n' ',' < ", outputPrefix, "_", phenoCol, "_colnames_subset_temp | sed 's/,$//') > ", outputPrefix, "_", phenoCol, "_subcols_temp")
      system(cmdb)

      phenoFiletemp <- paste0(outputPrefix, "_", phenoCol, "_subcols_temp")

      data <- data.table::fread(phenoFiletemp,
        header = T,
        stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F
      )
      # data = data.frame(ydat)

      file.remove(paste0(outputPrefix, "_", phenoCol, "_colnames_subset_temp"))
      file.remove(paste0(outputPrefix, "_", phenoCol, "_lineNum_temp"))
      file.remove(paste0(outputPrefix, "_", phenoCol, "_subcols_temp"))
    } else { # !isphenoFileLarge

      if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
        data <- data.table::fread(
          cmd = paste0(
            "gunzip -c ",
            phenoFile
          ), header = T, stringsAsFactors = FALSE,
          colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      } else {
        data <- data.table::fread(phenoFile,
          header = T,
          stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table = F, select = checkColList
        )
      }
    }

    file.remove(paste0(outputPrefix, "_", phenoCol, "_size_temp"))


    if (isRemoveZerosinPheno) {
      data <- data[which(data[, which(colnames(data) == phenoCol)] > 0), ]
      cat("Removing all zeros in the phenotype\n")
      if (nrow(data) == 0) {
        stop("ERROR: no samples are left after removing zeros in the phenotype\n")
      }
    }



    if (SampleIDIncludeFile != "") {
      if (!file.exists(SampleIDIncludeFile)) {
        stop("ERROR! SampleIDIncludeFile ", SampleIDIncludeFile, " does not exsit\n")
      } else {
        sampleIDInclude <- data.table::fread(SampleIDIncludeFile, header = F, stringsAsFactors = FALSE, colClasses = c("character"), data.table = F)
        sampleIDInclude <- as.vector(sampleIDInclude[!duplicated(sampleIDInclude), ])
        cat(length(sampleIDInclude), " non-duplicated sample IDs were found in SampleIDIncludeFile\n")
        data <- data[which(as.vector(data[, which(colnames(data) == sampleIDColinphenoFile)]) %in% sampleIDInclude), , drop = F]
        cat(nrow(data), " samples in sampleIDInclude have non-missing phenotypes and covariates\n")
      }
    }


    if (length(qCovarCol) > 0) {
      cat(qCovarCol, "are categorical covariates\n")
      if (!all(qCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in qCovarCol must be in covarColList\n")
      } else {
        for (q in qCovarCol) {
          data[, q] <- as.factor(data[, q])
        }
      }
    }

    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      if (!all(eCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in eCovarCol must be in covarColList\n")
      }
    }


    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      if (!all(sampleCovarCol %in% covarColList)) {
        stop("ERROR! all covariates in sampleCovarCol must be in covarColList\n")
      }
    }

    if (FemaleOnly | MaleOnly) {
      if (!sexCol %in% colnames(data)) {
        stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
      } else {
        if (FemaleOnly) {
          data <- data[which(data[, which(colnames(data) ==
            sexCol)] == FemaleCode), ]
          if (nrow(data) == 0) {
            stop(
              "ERROR! no samples in the phenotype are coded as ",
              FemaleCode, " in the column ", sexCol,
              "\n"
            )
          }
        } else if (MaleOnly) {
          data <- data[which(data[, which(colnames(data) ==
            sexCol)] == MaleCode), ]
          if (nrow(data) == 0) {
            stop(
              "ERROR! no samples in the phenotype are coded as ",
              MaleCode, " in the column ", sexCol, "\n"
            )
          }
        }
      }
    }


    print("HERERE2")

    if (length(covarColList) > 0) {
      formula <- paste0(phenoCol, "~", paste0(covarColList,
        collapse = "+"
      ))
      hasCovariate <- TRUE
    } else {
      formula <- paste0(phenoCol, "~ 1")
      hasCovariate <- FALSE
    }

    cat("formula is ", formula, "\n")
    formula.null <- as.formula(formula)
    mmat <- model.matrix(formula.null, data, na.action = NULL)
    mmat <- cbind(mmat, data[, which(colnames(data) == phenoCol), drop = F])
    colnames(mmat)[ncol(mmat)] <- phenoCol

    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")
      # check which sample-level covariates are categorical and record the names after factorizing in the data frame
      if (length(qCovarCol) > 0) {
        if (any(sampleCovarCol %in% qCovarCol)) {
          sampleCovarCol_q <- sampleCovarCol[which(sampleCovarCol %in% qCovarCol)]
          formula_sq <- paste0("~", paste0(sampleCovarCol_q, collapse = "+"))
          formula_sq.null <- as.formula(formula_sq)
          mmat_sq <- model.matrix(formula_sq.null, data, na.action = NULL)
          sampleCovarCol_q_names <- colnames(mmat_sq)[-1]
          rm(mmat_sq)
        } else {
          sampleCovarCol_q_names <- NULL
        }
      } else {
        sampleCovarCol_q_names <- NULL
      }
    }



    coln <- 1
    if (length(offsetCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == offsetCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- offsetCol
      coln <- coln + 1
    }

    if (length(varWeightsCol) > 0) {
      mmat <- cbind(mmat, data[, which(colnames(data) == varWeightsCol), drop = F])
      colnames(mmat)[ncol(mmat)] <- varWeightsCol

      coln <- coln + 1
    }



    if (length(covarColList) > 0) {
      if (length(qCovarCol) > 0) {
        covarColList <- colnames(mmat)[2:(ncol(mmat) - coln)]
        formula <- paste0(phenoCol, "~", paste0(covarColList, collapse = "+"))
        formula.null <- as.formula(formula)
      }
    }

    mmat$IID <- data[, which(sampleIDColinphenoFile == colnames(data))]
    if (cellIDColinphenoFile != "") {
      mmat$barcode <- data[, which(cellIDColinphenoFile == colnames(data))]
    }
    if (longlCol != "") {
      mmat$longlVar <- data[, which(longlCol == colnames(data))]
    }

    mmat_nomissing <- mmat[complete.cases(mmat), ]
    mmat_nomissing$IndexPheno <- seq(1, nrow(mmat_nomissing),
      by = 1
    )
    cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")

    if (length(varWeightsCol) > 0) {
      varWeights <- mmat_nomissing[, which(colnames(mmat_nomissing) == varWeightsCol)]
    } else {
      varWeights <- NULL
    }
    if (sparseGRMSampleIDFile != "") {
      # if((useSparseGRMtoFitNULL & !skipVarianceRatioEstimation) | useSparseGRMforVarRatio){
      sampleListwithGenov0 <- data.table::fread(sparseGRMSampleIDFile,
        header = F, , colClasses = c("character"), data.table = F
      )
      colnames(sampleListwithGenov0) <- c("IIDgeno")
      cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
      mmat_nomissing <- mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
      cat(nrow(mmat_nomissing), " samples who have non-missing phenotypes are also in the sparse GRM\n")
    }


    if (longlCol == "") {
      if (any(duplicated(mmat_nomissing$IID))) {
        cat("Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements\n")
      }
    } else {
      cat("Longitudinal variable ", longlCol, " is specified\n")
      if (!any(duplicated(mmat_nomissing$IID))) {
        stop("No duplicated sample IDs are detected in the phenotype file\n")
      }
    }


    if (!is.null(sampleListwithGeno)) {
      dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
        by.x = "IID", by.y = "IIDgeno"
      )
      dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]
      # dataMerge_sort = dataMerge[with(dataMerge, order(IndexPheno)),]
    } else {
      dataMerge_sort <- mmat_nomissing
      dataMerge_sort$IIDgeno <- dataMerge_sort$IID
    }

    print("Test")
    print(head(dataMerge_sort))
    # dataMerge_sort = dataMerge

    rm(mmat)
    rm(mmat_nomissing)
    gc()
    isSparseGRMIdentity <- FALSE
    if (useGRMtoFitNULL) {
      indicatorGenoSamplesWithPheno <- (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

      if (length(unique(dataMerge_sort$IIDgeno)) < length(unique(sampleListwithGeno$IIDgeno))) {
        cat(
          length(unique(sampleListwithGeno$IIDgeno)) - length(unique(dataMerge_sort$IIDgeno)),
          " samples in geno file do not have phenotypes\n"
        )
      }
      cat(length(unique(dataMerge_sort$IIDgeno)), " samples will be used for analysis\n")
    } else {
      indicatorGenoSamplesWithPheno <- rep(TRUE, nrow(dataMerge_sort))
    }

    if (any(duplicated(dataMerge_sort$IID))) {
      cat(nrow(dataMerge_sort), " observations will be used for analysis\n")
      # if(longlCol != ""){
      # dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(IndexGeno, longlVar)),]
      #  dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(longlVar)),]
      # }
      set_I_mat_inR(dataMerge_sort$IID)
      if (longlCol != "") {
        set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
      }
    } else {
      if (!useGRMtoFitNULL) {
        # stop("No duplicated IDs are observed in the phenotype file, so GRM must be used to fit the null model. Please set useGRMtoFitNULL=TRUE\n")
        cat("No duplicated IDs are observed in the phenotype file, so the identity matrix will be used as a sparse GRM will be used to fit the null model\n")
        isSparseGRMIdentity <- TRUE
        useSparseGRMtoFitNULL <- TRUE
        useGRMtoFitNULL <- TRUE
      }
    }
    set_useGRMtoFitNULL(useGRMtoFitNULL)
  }



  print("Test3")
  print(head(dataMerge_sort))

  if (traitType == "quantitative" & invNormalize) {
    cat(
      "Perform the inverse nomalization for ", phenoCol,
      "\n"
    )
    invPheno <- qnorm((rank(dataMerge_sort[, which(colnames(dataMerge_sort) ==
      phenoCol)], na.last = "keep") - 0.5) / sum(!is.na(dataMerge_sort[
      ,
      which(colnames(dataMerge_sort) == phenoCol)
    ])))
    dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] <- invPheno
  }
  print("Test4")
  print(head(dataMerge_sort))
  if (traitType == "binary" & (length(covarColList) > 0)) {
    out_checksep <- checkPerfectSep(formula.null,
      data = dataMerge_sort,
      minCovariateCount
    )
    covarColList <- covarColList[!(covarColList %in% out_checksep)]
    formula <- paste0(phenoCol, "~", paste0(covarColList,
      collapse = "+"
    ))
    formula.null <- as.formula(formula)
    if (length(covarColList) == 1) {
      hasCovariate <- FALSE
    } else {
      hasCovariate <- TRUE
    }
    dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in%
      out_checksep)]
  }
  if (!hasCovariate) {
    print("No covariate is includes so isCovariateOffset = FALSE")
    isCovariateOffset <- FALSE
  }

  if (isCovariateTransform & hasCovariate) {
    cat("qr transformation has been performed on covariates\n")
    out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
    formulaNewList <- c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
    if (length(out.transform$Param.transform$X_name) > 1) {
      for (i in c(2:length(out.transform$Param.transform$X_name))) {
        formulaNewList <- c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
      }
    }
    formulaNewList <- paste0(formulaNewList, collapse = "")
    formulaNewList <- paste0(formulaNewList, "-1")
    formula.new <- as.formula(paste0(formulaNewList, collapse = ""))
    data.new <- as.data.frame(cbind(out.transform$Y, out.transform$X1))
    colnames(data.new) <- c(phenoCol, out.transform$Param.transform$X_name)
    cat("colnames(data.new) is ", colnames(data.new), "\n")
    cat(
      "out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr),
      "\n"
    )

    if (length(offsetCol) > 0) {
      data.new <- cbind(data.new, dataMerge_sort[which(colnames(dataMerge_sort) == offsetCol), ])
      colnames(data.new)[ncol(data.new)] == offsetCol
    }
  } else {
    formula.new <- formula.null
    data.new <- dataMerge_sort
    out.transform <- NULL
  }

  if (traitType == "binary") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = binomial, weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = binomial, weights = varWeights
      )
    }
  } else if (traitType == "quantitative") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = gaussian(link = "identity"), weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = gaussian(link = "identity"), weights = varWeights
      )
    }
  } else if (traitType == "count") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = "poisson", weights = varWeights
      )
    } else {
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = "poisson", weights = varWeights
      )
    }
  } else if (traitType == "count_nb") {
    if (length(offsetCol) == 0) {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = NegBin(), weights = varWeights
      )
    } else {
      print(head(data.new))
      offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
      modwitcov <- glm(formula.new,
        offset = offsetColVal, data = data.new,
        family = NegBin(), weights = varWeights
      )
    }
  }
  mmat <- model.matrix(formula.new, data = data.new, na.action = NULL)

  if (isCovariateOffset) {
    covoffset <- mmat[, -1, drop = F] %*% modwitcov$coefficients[-1]
    print("isCovariateOffset=TRUE, so fixed effects coefficnets won't be estimated.")
    formula.new.withCov <- formula.new
    formula_nocov <- paste0(phenoCol, "~ 1")
    formula.new <- as.formula(formula_nocov)
    hasCovariate <- FALSE
  } else {
    covoffset <- rep(0, nrow(data.new))
  }

  data.new$covoffset <- covoffset


  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!isSparseGRMIdentity) {
      getsubGRM_orig(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID)
    } else {
      print(length(dataMerge_sort$IIDgeno))
      sparseGRM <- Matrix:::sparseMatrix(i = as.vector(1:nrow(data.new)), j = as.vector(1:nrow(data.new)), x = rep(1, nrow(data.new)), symmetric = TRUE)
      setupSparseGRM_new(sparseGRM)
    }
    # getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID, dataMerge_sort$longlVar)
    # m4 = gen_sp_v2(sparseGRMtest)
    # cat("Setting up sparse GRM using ", sparseGRMFile, " and ", sparseGRMSampleIDFile, "\n")
    # cat("Dimension of the sparse GRM is ", dim(m4), "\n")
    # A = summary(m4)
    # locationMatinR = rbind(A$i - 1, A$j - 1)
    # valueVecinR = A$x
    # setupSparseGRM(dim(m4)[1], locationMatinR, valueVecinR)
    # setupSparseGRM_new(sparseGRMtest)
    # rm(sparseGRMtest)
    gc()
  }

  # allow for multiple variance components
  # set_Vmat_vec(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID, dataMerge_sort$longlVar)
  set_Vmat_vec_orig(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID)

  numofV <- get_numofV()

  print(dataMerge_sort$IID[1:200])
  print(any(duplicated(dataMerge_sort$IID)))

  cat("numofV ", numofV, "\n")
  if (any(duplicated(dataMerge_sort$IID))) {
    print("HERE")
    if (longlCol == "") {
      print("HERE1")
      print(useGRMtoFitNULL)
      if (useGRMtoFitNULL) {
        print("HERE2")
        num_Kmat <- numofV + 3
        cat("num_Kmat ", num_Kmat, "\n")
      } else {
        num_Kmat <- numofV + 2
      }
      # k = num_Kmat + 3
    } else {
      if (useGRMtoFitNULL) {
        num_Kmat <- 7 + numofV * 3
      } else {
        num_Kmat <- 4 + numofV * 3
      }
      # k = num_Kmat + 2
    }
  } else {
    if (useGRMtoFitNULL) {
      num_Kmat <- numofV + 2
    } else {
      num_Kmat <- numofV + 1
    }
    # k = 2
  }

  k <- num_Kmat


  set_num_Kmat(num_Kmat)
  cat("num_Kmat ", num_Kmat, "\n")


  if (longlCol != "") {
    covarianceIdxMat <- set_covarianceidx_Mat()
  } else {
    covarianceIdxMat <- NULL
  }

  if (!skipVarianceRatioEstimation) {
    isVarianceRatioinGeno <- TRUE
    if (isCateVarianceRatio) {
      minMAC_varRatio <- min(cateVarRatioMinMACVecExclude)
      maxMAC_varRatio <- max(cateVarRatioMaxMACVecInclude)
      cat("Categorical variance ratios will be estimated. Please make sure there are at least 200 markers in each MAC category.\n")
    } else {
      minMAC_varRatio <- 20
      maxMAC_varRatio <- -1 # will randomly select markers from the plink file and leave them out when constructing GRM
    }
    setminMAC_VarianceRatio(minMAC_varRatio, maxMAC_varRatio, isVarianceRatioinGeno)
  }

  # set up parameters
  if (minMAFforGRM > 0) {
    cat(
      "Markers in the Plink file with MAF < ", minMAFforGRM,
      " will be removed before constructing GRM\n"
    )
  }
  if (maxMissingRateforGRM > 0) {
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  setminMAFforGRM(minMAFforGRM)
  setmaxMissingRateforGRM(maxMissingRateforGRM)



  if (traitType == "binary") {
    cat(phenoCol, " is a binary trait\n")
    uniqPheno <- sort(unique(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)]))
    if (uniqPheno[1] != 0 | uniqPheno[2] != 1) {
      stop("ERROR! phenotype value needs to be 0 or 1 \n")
    }
    print("formula.new")
    print(formula.new)
    print("head(data.new)")
    print(head(data.new))
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = binomial, weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = binomial, weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = binomial, weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = binomial, weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = binomial, weights = varWeights)
      }
    }
  } else if (traitType == "quantitative") {
    cat(phenoCol, " is a quantitative trait\n")
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = gaussian(link = "identity"), weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = gaussian(link = "identity"), weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = gaussian(link = "identity"), weights = varWeights)
      }
    }
  } else if (traitType == "count") {
    cat(phenoCol, " is a count trait\n")
    # print("before remove zeros")
    # print(dim(dataMerge_sort))
    # if(isRemoveZerosinPheno){
    #    dataMerge_sort = dataMerge_sort[which(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] > 0),]
    #    cat("Removing all zeros in the phenotype\n")
    #    if(nrow(dataMerge_sort) == 0){
    #        stop("ERROR: no samples are left after removing zeros in the phenotype\n")

    #    }
    # }
    # print("after remove zeros")
    # print(dim(dataMerge_sort))
    miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
    if (miny < 0) {
      stop("ERROR! phenotype value needs to be non-negative \n")
    }

    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = "poisson", weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = "poisson", weights = varWeights)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = "poisson", weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = "poisson", weights = varWeights
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = "poisson", weights = varWeights)
      }
    }
  } else if (traitType == "count_nb") {
    cat(phenoCol, " is a count_nb trait\n")
    if (isRemoveZerosinPheno) {
      dataMerge_sort <- dataMerge_sort[which(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] > 0), ]
      cat("Removing all zeros in the phenotype\n")
      if (nrow(dataMerge_sort) == 0) {
        stop("ERROR: no samples are left after removing zeros in the phenotype\n")
      }
    }
    miny <- min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
    if (miny < 0) {
      stop("ERROR! phenotype value needs to be non-negative \n")
    }
    if (!isCovariateOffset) {
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new, data = data.new, family = NegBin(), weights = varWeights)
      } else {
        offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = NegBin(), weights = varWeights)
      }

      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = NegBin(), weights = varWeights)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      if (length(offsetCol) == 0) {
        fit0 <- glm(formula.new,
          data = data.new, offset = covoffset,
          family = NegBin()
        )
      } else {
        offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
        fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = NegBin(), weights = varWeights)
      }
    }
  }


  cat("glm:\n")
  print(fit0)
  obj.noK <- NULL


  # if(length(fit0$y) > 200000){
  #  isStoreSigma = FALSE
  # }else{
  #  isStoreSigma = TRUE
  # }
  print("isStoreSigma")
  print(isStoreSigma)
  # set_store_sigma(isStoreSigma)

  if (!skipModelFitting) {
    # setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
    cat("Start fitting the NULL GLMM\n")
    t_begin <- proc.time()
    print(t_begin)


    tau <- rep(0, k)
    fixtau <- rep(0, k)
    # tauInit = tau

    set_isSparseGRM(useSparseGRMtoFitNULL)
    set_useGRMtoFitNULL(useGRMtoFitNULL)

    if (traitType != "count_nb") {
      system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = tau, fixtau = fixtau, maxiter = maxiter,
        tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL, isSparseGRMIdentity = isSparseGRMIdentity
      ))
      modglmm$obj.glm.null$model <- as.data.frame(modglmm$obj.glm.null$model)
    } else {
      system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV_NB(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = tau, fixtau = fixtau, maxiter = maxiter,
        tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL, isSparseGRMIdentity = isSparseGRMIdentity
      ))

      data.new$y <- modglmm$y
      varWeights <- modglmm$varWeights
      if (!isCovariateOffset) {
        if (length(offsetCol) == 0) {
          fit0 <- glm(formula.new, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
        } else {
          offsetColVal <- data.new[, which(colnames(data.new) == offsetCol)]
          fit0 <- glm(formula.new, data = data.new, offset = offsetColVal, family = gaussian(link = "identity"), weights = varWeights)
        }
        Xorig <- NULL
      } else {
        fit0orig <- glm(formula.new.withCov, data = data.new, family = gaussian(link = "identity"), weights = varWeights)
        Xorig <- model.matrix(fit0orig)
        rm(fit0orig)
        gc()
        if (length(offsetCol) == 0) {
          fit0 <- glm(formula.new,
            data = data.new, offset = covoffset,
            family = gaussian(link = "identity"), weights = varWeights
          )
        } else {
          offsetTotal <- covoffset + data.new[, which(colnames(data.new) == offsetCol)]
          fit0 <- glm(formula.new, data = data.new, offset = offsetTotal, family = gaussian(link = "identity"), weights = varWeights)
        }
      }

      modglmm$obj.glm.null <- fit0
    }


    # spSigma_final = getSparseSigma_new()
    # modglmm$spSigma = spSigma_final
    # rm(spSigma_final)
    if (traitType != "count_nb") {
      for (x in names(modglmm$obj.glm.null)) {
        attr(modglmm$obj.glm.null[[x]], ".Environment") <- c()
      }
    }
    # modglmm$offset = covoffset
    if (isCovariateOffset) {
      modglmm$offset <- covoffset
    } else {
      if (hasCovariate) {
        data.new.X <- model.matrix(fit0)[, -1, drop = F]
        print(head(data.new))
        print(head(data.new.X))
        print(head(modglmm$coefficients[-1]))
        modglmm$offset <- data.new.X %*% (as.vector(modglmm$coefficients[-1]))
        if (LOCO) {
          for (j in 1:22) {
            if (modglmm$LOCOResult[[j]]$isLOCO) {
              modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
            }
          }
        }
      } else {
        modglmm$offset <- covoffset
      }
    }

    if (length(eCovarCol) > 0) {
      cat(eCovarCol, "are environmental covariates\n")
      modglmm$eMat <- data.new[, which(colnames(data.new) %in% eCovarCol), drop = F]
      for (em in 1:ncol(modglmm$eMat)) {
        modglmm$eMat[, em] <- (modglmm$eMat[, em] - mean(modglmm$eMat[, em])) / (sd(modglmm$eMat[, em]))
      }
    }

    if (length(sampleCovarCol) > 0) {
      cat(sampleCovarCol, "are sample-level covariates\n")

      sampleCovarCol <- c(sampleCovarCol, sampleCovarCol_q_names)
      modglmm$sampleXMat <- modglmm$X[, which(colnames(modglmm$X) %in% sampleCovarCol), drop = F]
      modglmm$sampleXMat <- cbind(modglmm$X[, 1], modglmm$sampleXMat)
      uniqsampleind <- which(!duplicated(modglmm$sampleID))
      modglmm$sampleXMat <- modglmm$sampleXMat[uniqsampleind, ]
    }


    # if((skipVarianceRatioEstimation & useSparseGRMtoFitNULL)){
    family <- fit0$family
    eta <- modglmm$linear.predictors
    mu <- modglmm$fitted.values
    mu.eta <- family$mu.eta(eta)
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2
    W <- W * modglmm$varWeights
    tauVecNew <- modglmm$theta
    Sigma_iX <- getSigma_X_multiV(W, tauVecNew, modglmm$X, maxiterPCG, tolPCG, LOCO = FALSE)
    if (!isShrinkModelOutput) {
      Sigma_iXXSigma_iX <- Sigma_iX %*% solve(crossprod(modglmm$X, Sigma_iX))
      modglmm$Sigma_iXXSigma_iX <- Sigma_iXXSigma_iX
    }
    # }

    modglmm$useSparseGRMforVarRatio <- useSparseGRMforVarRatio


    # if(any(duplicated(modglmm$sampleID))){
    if (useGRMtoFitNULL) {
      modglmm$tauVal_sp <- modglmm$theta[3]
    } else {
      modglmm$tauVal_sp <- modglmm$theta[2]
    }
    # }

    # if(FALSE){
    # if(length(fit0$y) <= 10000){

    cat("isStoreSigma is ", isStoreSigma, "\n")
    if (isStoreSigma) {
      # family = fit0$family
      # eta = modglmm$linear.predictors
      # mu = modglmm$fitted.values
      # mu.eta = family$mu.eta(eta)
      # sqrtW = mu.eta/sqrt(family$variance(mu))
      # W = sqrtW^2
      # W = W * modglmm$varWeights;
      # tauVecNew = modglmm$theta
      # Matrix::writeMM(spSigma, file="/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/saigeqtl_manuscript/realdata/step1/LSM10.CD4_NC.spSigma.mtx")
      # if(any(duplicated(dataMerge_sort$IID))){
      modglmm$spSigma <- gettI_Sigma_I_multiV(W, tauVecNew, maxiterPCG, tolPCG, LOCO = FALSE)
      # }else{
      # 	gen_sp_Sigma_multiV(W, tauVecNew)
      # 	spSigma = get_sp_Sigma_to_R()
      # 	SigmaMat_sp = chol2inv(chol(spSigma))
      # 	modglmm$spSigma = SigmaMat_sp
      # }
    }
    # }
    # save(modglmm, file = modelOut)
    tau <- modglmm$theta
    alpha0 <- modglmm$coefficients


    if (!is.null(out.transform) & is.null(fit0$offset)) {
      coef.alpha <- Covariate_Transform_Back(alpha0, out.transform$Param.transform)
      modglmm$coefficients <- coef.alpha
    }


    if (LOCO & isLowMemLOCO) {
      modglmm$LOCOResult <- NULL
      modglmm$LOCO <- FALSE
      chromosomeStartIndexVec <- modglmm$chromosomeStartIndexVec
      chromosomeEndIndexVec <- modglmm$chromosomeEndIndexVec
      modglmm$chromosomeStartIndexVec <- NULL
      modglmm$chromosomeEndIndexVec <- NULL
      modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
      fastSave(modglmm, file = modelOut)
      modglmm$LOCO <- TRUE
      modglmm$Y <- NULL
      eta0 <- modglmm$linear.predictors
      modglmm$linear.predictors <- NULL
      modglmm$coefficients <- NULL
      modglmm$cov <- NULL
      modglmm$fitted.values <- NULL
      modglmm$residuals <- NULL
      modglmm$obj.noK <- NULL
      offset0 <- modglmm$offset
      modglmm$offset <- NULL
      y <- fit0$y
      gc()
      # save(modglmm, file = modelOut)
      set_Diagof_StdGeno_LOCO()
      modglmm$LOCOResult <- list()
      for (j in 1:22) {
        startIndex <- chromosomeStartIndexVec[j]
        endIndex <- chromosomeEndIndexVec[j]
        if (!is.na(startIndex) && !is.na(endIndex)) {
          cat("leave chromosome ", j, " out\n")
          setStartEndIndex(startIndex, endIndex, j - 1)

          re.coef_LOCO <- Get_Coef_multiV(y, X = model.matrix(fit0), tau, family = fit0$family, alpha = alpha0, eta = eta0, offset = offset0, verbose = TRUE, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = TRUE, var_weights = var_weights)
          cov <- re.coef_LOCO$cov
          alpha <- re.coef_LOCO$alpha
          eta <- re.coef_LOCO$eta
          Y <- re.coef_LOCO$Y
          mu <- re.coef_LOCO$mu
          if (family$family == "binomial") {
            mu2 <- mu * (1 - mu)
          } else if (family$family == "poisson") {
            mu2 <- mu
          } else if (family$family == "gaussian") {
            mu2 <- rep((1 / (tau[1])), length(res))
          } else if (traitType == "count_nb") {
            # mu2 = fit0$family$variance(mu)
            mu2 <- rep((1 / (tau[1])), length(res))
          }
          res <- y - mu

          if (!is.null(out.transform) & is.null(fit0$offset)) {
            coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
          } else {
            coef.alpha <- alpha
          }

          mu2_rescaled <- mu2 * var_weights
          mu_rescaled <- mu * var_weights


          if (!isCovariateOffset) {
            obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
          } else {
            obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
          }
          modglmm$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)
          if (!isCovariateOffset & hasCovariate) {
            data.new.X <- model.matrix(fit0)[, -1, drop = F]
            modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
          }
          modelOutbychr <- paste(c(outputPrefix, "_chr", j, ".rda"), collapse = "")
          if (j != 22) {
            for (j1 in (j + 1):22) {
              modglmm$LOCOResult[[j1]] <- list(NULL)
            }
          }
          fastSave(modglmm, file = modelOutbychr)
          modglmm$LOCOResult[[j]] <- list(NULL)
          gc()
        } else {
          modglmm$LOCOResult[[j]] <- list(isLOCO = FALSE)
        }
      }
      gc()
      # modelOut_nonauto = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")
    } else {
      # b = as.numeric(factor(dataMerge_sort$IID, levels =  unique(dataMerge_sort$IID)))
      # I_mat = Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
      # I_mat = 1.0 * I_mat
      # modglmm$I_longl_mat = I_mat
      # modglmm$I_longl_vec = b - 1
      # modglmm$T_longl_mat = I_mat * (dataMerge_sort$longlVar)
      modglmm$T_longl_vec <- dataMerge_sort$longlVar
    }


    t_end <- proc.time()
    print(t_end)
    cat("t_end - t_begin, fitting the NULL model took\n")
    print(t_end - t_begin)


    # if(bedFile != "" & !useGRMtoFitNULL){
    if (bedFile != "") {
      subSampleInGeno <- dataMerge_sort$IndexGeno
      if (is.null(dataMerge_sort$IndexGeno)) {
        subSampleInGeno <- dataMerge_sort$IndexPheno
      }

      print(subSampleInGeno[1:1000])
      print(head(dataMerge_sort))
      print("HEREHRE")

      subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]

      # setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
      setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    }
  } else {
    cat("Skip fitting the NULL GLMM\n")
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
    load(modelOut)
    if (is.null(modglmm$LOCO)) {
      modglmm$LOCO <- FALSE
    }

    # need check
    subSampleInGeno <- dataMerge_sort$IndexGeno
    if (is.null(dataMerge_sort$IndexGeno)) {
      subSampleInGeno <- dataMerge_sort$IndexPheno
    }

    print(subSampleInGeno[1:1000])
    print(head(dataMerge_sort))
    print("HEREHRE")

    subSampleInGeno_unique <- subSampleInGeno[!duplicated(subSampleInGeno)]

    # setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    setgeno(bedFile, bimFile, famFile, subSampleInGeno_unique, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)


    # setgeno(bedFile, bimFile, famFile, dataMerge_sort$IndexGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    tau <- modglmm$theta

    if (any(duplicated(modglmm$sampleID))) {
      set_I_mat_inR(modglmm$sampleID)
      # if(longlCol != ""){
      #        set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
      # }
    }
    set_dup_sample_index(as.numeric(factor(modglmm$sampleID, levels = unique(modglmm$sampleID))))
    # setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
  }

  if (!skipVarianceRatioEstimation) {
    if (LOCO) {
      MsubIndVec <- getQCdMarkerIndex()
      print(length(MsubIndVec))
      chrVec <- data.table::fread(bimFile, header = F)[, 1]
      print(length(chrVec))
      chrVec <- chrVec[which(MsubIndVec == TRUE)]
      updatechrList <- updateChrStartEndIndexVec(chrVec)
      LOCO <- updatechrList$LOCO
      chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
      chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
      set_Diagof_StdGeno_LOCO()
      load(modelOut)
    }
    cat("Start estimating variance ratios\n")
    extractVarianceRatio_multiV(
      obj.glmm.null = modglmm,
      obj.glm.null = fit0, maxiterPCG = maxiterPCG,
      tolPCG = tolPCG, numMarkers = numMarkersForVarRatio, varRatioOutFile = varRatioFile,
      ratioCVcutoff = ratioCVcutoff, testOut = SPAGMMATOut,
      bedFile = bedFile, bimFile = bimFile, famFile = famFile, chromosomeStartIndexVec = chromosomeStartIndexVec,
      chromosomeEndIndexVec = chromosomeEndIndexVec,
      isCateVarianceRatio = isCateVarianceRatio, cateVarRatioIndexVec = cateVarRatioIndexVec,
      useSparseGRMforVarRatio = useSparseGRMforVarRatio, sparseGRMFile = sparseGRMFile,
      sparseGRMSampleIDFile = sparseGRMSampleIDFile,
      numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
      relatednessCutoff = relatednessCutoff, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL,
      nThreads = nThreads, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      minMAFforGRM = minMAFforGRM, isDiagofKinSetAsOne = isDiagofKinSetAsOne,
      includeNonautoMarkersforVarRatio = includeNonautoMarkersforVarRatio, isStoreSigma = isStoreSigma, useGRMtoFitNULL = useGRMtoFitNULL
    )
  } else {
    cat("Skip estimating variance ratios\n")
  }
  closeGenoFile_plink()

  # clean up saved model (as in ReadModel)
  if (isShrinkModelOutput) {
    modglmm$Y <- NULL
    modglmm$linear.predictors <- NULL
    modglmm$coefficients <- NULL
    modglmm$cov <- NULL
    if (sum(duplicated(modglmm$sampleID)) > 0) {
      modglmm$obj.noK$Sigma_iXXSigma_iX <- matrix(1)
      modglmm$X <- NULL
      if (is.null(modglmm$eMat)) {
        modglmm$obj.noK$XV <- NULL
        modglmm$obj.noK$XVX <- NULL
        modglmm$obj.noK$XXVX_inv <- NULL
        modglmm$obj.noK$XVX_inv <- NULL
        modglmm$obj.noK$XVX_inv_XV <- NULL
      }
    }
  }
  fastSave(modglmm, file = modelOut)

  if (isExportResiduals) {
    b <- as.numeric(factor(modglmm$sampleID, levels = unique(modglmm$sampleID)))
    I_mat <- 1.0 * Matrix::sparseMatrix(i = seq_along(b), j = b, x = rep(1, length(b)))
    res_sample <- as.vector(t(I_mat) %*% modglmm$residuals)
    data.table::fwrite(
      data.frame(sampleID = unique(modglmm$sampleID), residuals = res_sample),
      paste0(outputPrefix, ".sample.residuals.txt"),
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      na = "NA"
    )
    data.table::fwrite(
      data.frame(barcode = modglmm$barcode, residuals = modglmm$residuals),
      paste0(outputPrefix, ".residuals.txt"),
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      na = "NA"
    )
  }
}



extractVarianceRatio_multiV <- function(obj.glmm.null,
                                        obj.glm.null,
                                        maxiterPCG = 500,
                                        tolPCG = 0.01,
                                        numMarkers,
                                        varRatioOutFile,
                                        ratioCVcutoff,
                                        testOut,
                                        bedFile,
                                        bimFile,
                                        famFile,
                                        chromosomeStartIndexVec,
                                        chromosomeEndIndexVec,
                                        isCateVarianceRatio,
                                        cateVarRatioIndexVec,
                                        useSparseGRMforVarRatio,
                                        sparseGRMFile,
                                        sparseGRMSampleIDFile,
                                        numRandomMarkerforSparseKin,
                                        relatednessCutoff,
                                        useSparseGRMtoFitNULL,
                                        nThreads,
                                        cateVarRatioMinMACVecExclude,
                                        cateVarRatioMaxMACVecInclude,
                                        minMAFforGRM,
                                        isDiagofKinSetAsOne,
                                        includeNonautoMarkersforVarRatio,
                                        isStoreSigma = FALSE,
                                        useGRMtoFitNULL = TRUE) {
  obj.noK <- obj.glmm.null$obj.noK
  if (file.exists(testOut)) {
    file.remove(testOut)
  }
  bimPlink <- data.table::fread(bimFile, header = F, data.table = FALSE)
  if (sum(sapply(bimPlink[, 1], is.numeric)) != nrow(bimPlink)) {
    stop("ERROR: chromosome column in plink bim file is no numeric!\n")
  }

  family <- obj.glm.null$family
  print(family)
  eta <- obj.glmm.null$linear.predictors
  mu <- obj.glmm.null$fitted.values
  mu.eta <- family$mu.eta(eta)

  # var_weights = weights(obj.glm.null)
  var_weights <- obj.glmm.null$varWeights
  # if(!is.null(var_weights)){
  sqrtW <- mu.eta / sqrt(family$variance(mu))
  # }else{
  #  sqrtW = mu.eta/sqrt(family$variance(mu))
  # }


  print("mu[1:20]")
  print(mu[1:20])
  W <- sqrtW^2 ## (mu*(1-mu) for binary)
  W <- W * var_weights
  print("mu[1:20]")
  print(W[1:20])
  tauVecNew <- obj.glmm.null$theta

  # isStoreSigma=FALSE
  # if(isStoreSigma){
  #       gen_sp_Sigma_multiV(W, tauVecNew)
  # }
  X <- obj.glmm.null$X


  set_isSparseGRM(useSparseGRMtoFitNULL)
  set_useGRMtoFitNULL(useGRMtoFitNULL)

  # if(!useGRMtoFitNULL){
  # if(useSparseGRMforVarRatio){
  # }
  # useSparseGRMforVarRatio = FALSE
  # }

  Sigma_iX_noLOCO <- getSigma_X_multiV(W, tauVecNew, X, maxiterPCG, tolPCG, LOCO = FALSE)


  y <- obj.glmm.null$y

  if (any(duplicated(obj.glmm.null$sampleID))) {
    dupSampleIndex <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  }


  ## randomize the marker orders to be tested
  if (FALSE) {
    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
      sparseSigma <- getSparseSigma(
        bedFile = bedFile, bimFile = bimFile, famFile = famFile,
        outputPrefix = varRatioOutFile,
        sparseGRMFile = sparseGRMFile,
        sparseGRMSampleIDFile = sparseGRMSampleIDFile,
        numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
        relatednessCutoff = relatednessCutoff,
        minMAFforGRM = minMAFforGRM,
        nThreads = nThreads,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        obj.glmm.null = obj.glmm.null,
        W = W, tauVecNew = tauVecNew
      )
      if (length(tauVecNew) > 2) {
        sparseSigma <- sparseSigma + getProdTauKmat(tauVecNew[3:length(tauVecNew)])
      }
    }
  }


  mMarkers <- gettotalMarker()
  listOfMarkersForVarRatio <- list()
  MACvector <- getMACVec()
  isVarianceRatioinGeno <- getIsVarRatioGeno()

  if (isVarianceRatioinGeno) {
    MACvector_forVarRatio <- getMACVec_forVarRatio()
    Indexvector_forVarRatio <- getIndexVec_forVarRatio()
    cat("length(MACvector): ", length(MACvector), "\n")
    cat("length(MACvector_forVarRatio): ", length(MACvector_forVarRatio), "\n")

    if (length(MACvector_forVarRatio) > 0) {
      # MACdata = data.frame(MACvector = c(MACvector, MACvector_forVarRatio), geno_ind = c(rep(0, length(MACvector)), rep(1, length(MACvector_forVarRatio))), indexInGeno = c(seq(1,length(MACvector)), seq(1,length(MACvector_forVarRatio))))
      MACdata <- data.frame(MACvector = MACvector_forVarRatio, geno_ind = rep(1, length(MACvector_forVarRatio)), indexInGeno = seq(1, length(MACvector_forVarRatio)))
    } else {
      stop("No markers were found for variance ratio estimation. Please make sure there are at least 200 markers in each MAC category\n")
    }
  } else {
    MACdata <- data.frame(MACvector = MACvector, geno_ind = rep(0, length(MACvector)), indexInGeno = seq(1, length(MACvector)))
  }

  if (!isCateVarianceRatio) {
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex <- 1:nrow(MACdata)
    listOfMarkersForVarRatio[[1]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec <- c(1)
  } else {
    cat("Categorical variance ratios will be estimated.\n")

    if (is.null(cateVarRatioIndexVec)) {
      cateVarRatioIndexVec <- rep(1, length(cateVarRatioMinMACVecExclude))
    }
    numCate <- length(cateVarRatioIndexVec)
    for (i in 1:(numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[i] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[i])
      listOfMarkersForVarRatio[[i]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    }

    if (length(cateVarRatioMaxMACVecInclude) == (numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate])
    } else {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }

    listOfMarkersForVarRatio[[numCate]] <- sample(MACindex, size = length(MACindex), replace = FALSE)

    for (k in 1:length(cateVarRatioIndexVec)) {
      if (k <= length(cateVarRatioIndexVec) - 1) {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], "\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      } else {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }
  } # if(!isCateVarianceRatio){



  b <- as.numeric(factor(obj.glmm.null$sampleID, levels = unique(obj.glmm.null$sampleID)))
  I_mat <- Matrix::sparseMatrix(i = seq_along(b), j = b, x = rep(1, length(b)))

  freqVec <- getAlleleFreqVec()
  Nnomissing <- length(mu)
  varRatioTable <- NULL


  # uniqsampleind = which(!duplicated(obj.glmm.null$sampleID))
  # var_weights_sample = var_weights[uniqsampleind]


  Vsample0 <- as.vector(crossprod(obj.noK$V, I_mat))
  Xsample0 <- obj.glmm.null$sampleXMat
  XVsample0 <- t(Xsample0 * Vsample0)
  XVXsample0 <- crossprod(Xsample0, t(XVsample0))
  XVXsample_inv0 <- solve(XVXsample0)
  XXVXsample_inv0 <- Xsample0 %*% XVXsample_inv0
  XVX_inv_XVsample0 <- XXVXsample_inv0 * Vsample0
  # XVsample0_e_list = list()
  # for(ne in 1:ncol(obj.glmm.null$eMat)){
  #  evec = obj.glmm.null$eMat[,ne]
  #  Vsample0_e = as.vector(t(obj.noK$V * evec) %*% I_mat)
  #
  #    XVsample0_e =  t(Xsample0 * (Vsample0_e))
  #    XVsample0_e_list[[ne]] = XVsample0_e
  #  }

  for (k in 1:length(listOfMarkersForVarRatio)) {
    if (cateVarRatioIndexVec[k] == 1) {
      numMarkers0 <- numMarkers
      varRatio_sparseGRM_vec <- NULL
      varRatio_NULL_vec <- NULL
      varRatio_NULL_sample_vec <- NULL
      varRatio_NULL_noXadj_vec <- NULL


      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_mat <- NULL
        varRatio_NULL_eg_vec <- NULL
        varRatio_sparse_eg_mat <- NULL
        varRatio_sparse_eg_vec <- NULL
      }

      indexInMarkerList <- 1
      numTestedMarker <- 0
      ratioCV <- ratioCVcutoff + 0.1

      while (ratioCV > ratioCVcutoff) {
        while (numTestedMarker < numMarkers0) {
          macdata_i <- listOfMarkersForVarRatio[[k]][indexInMarkerList]
          i <- (MACdata$indexInGeno)[macdata_i]
          genoInd <- (MACdata$geno_ind)[macdata_i]
          cat(i, "th marker in geno ", genoInd, "\n")
          cat("MAC: ", (MACdata$MACvector)[macdata_i], "\n")
          if (genoInd == 0) {
            G0 <- Get_OneSNP_Geno(i - 1)
          } else if (genoInd == 1) {
            G0 <- Get_OneSNP_Geno_forVarRatio(i - 1)
          }

          # if(sum(duplicated(obj.glmm.null$sampleID)) > 0){
          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          G0sample <- G0
          print("length(G0)   aaaaa")
          print(length(G0))
          cat("G0", G0[1:10], "\n")
          print(dim(I_mat))
          print(length(G0sample))
          G0 <- as.numeric(I_mat %*% G0sample)
          # }
          # print("length(G0)   bbbb")
          # print(length(G0))
          cat("G0", G0[1:10], "\n")

          CHR <- bimPlink[Indexvector_forVarRatio[i] + 1, 1]
          cat("CHR ", CHR, "\n")
          print(bimPlink[Indexvector_forVarRatio[i] + 1, ])
          # if(sum(G0)/(2*Nnomissing) > 0.5){
          if (sum(G0) / (2 * length(G0)) > 0.5) {
            G0 <- 2 - G0
          }
          print("length(G0)")
          print(length(G0))
          # if(any(duplicated(obj.glmm.null$sampleID))){
          # 	G0 = G0[dupSampleIndex]
          # }
          NAset <- which(G0 == 0)
          AC <- sum(G0)
          print("length(NAset)")
          print(length(NAset))
          indexInMarkerList <- indexInMarkerList + 1
          if ((CHR >= 1 & CHR <= 22 & AC > 0 & AC < length(G0)) | includeNonautoMarkersforVarRatio) {
            AF <- AC / (2 * Nnomissing)
            if (CHR >= 1 & CHR <= 22) {
              autoMarker <- TRUE
            } else {
              autoMarker <- FALSE
            }

            G <- G0 - obj.noK$XXVX_inv %*% (obj.noK$XV %*% G0) # G1 is X adjusted
            # g = G/sqrt(AC)
            # q = innerProduct(g * sqrt(var_weights),y)
            # q = innerProduct(G,y)
            # eta = obj.glmm.null$linear.predictors
            # mu = obj.glmm.null$fitted.values
            # mu.eta = family$mu.eta(eta)
            # sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
            # W = sqrtW^2
            # W = W * var_weights
            # W = W * var_weights
            # print("W[1:10]")
            # print(W[1:10])
            # print("mu.eta[1:10]")
            # print(mu.eta[1:10])
            # print("obj.glm.null$family$variance(mu)[1:10]")
            # print(obj.glm.null$family$variance(mu)[1:10])
            # print("mu[1:100]")
            # print(mu[1:100])
            # print("y[1:100]")
            # print(y[1:100])


            set_isSparseGRM(useSparseGRMtoFitNULL)
            set_useGRMtoFitNULL(useGRMtoFitNULL)


            Sigma_iG <- getSigma_G_multiV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO = FALSE)
            Sigma_iX <- Sigma_iX_noLOCO

            var1 <- crossprod(G, Sigma_iG) - crossprod(G, Sigma_iX) %*% solve(crossprod(X, Sigma_iX)) %*% crossprod(X, Sigma_iG)
            cat("AC ", AC, "\n")
            S <- innerProduct(G, obj.glmm.null$residuals * var_weights)
            cat("S is ", S, "\n")
            cat("var1 is ", var1, "\n")
            p_exact <- pchisq(S^2 / var1, df = 1, lower.tail = F)
            cat("p_exact ", p_exact, "\n")
            # res_sample = as.vector(obj.glmm.null$residuals %*% I_mat)



            if (!is.null(obj.glmm.null$eMat)) {
              var1GE_vec <- NULL
              var2sparseGE_vec <- NULL
              getildeMat <- NULL
              getilde_sample0_Mat <- NULL
              for (ne in 1:ncol(obj.glmm.null$eMat)) {
                evec <- obj.glmm.null$eMat[, ne]

                # evec = rnorm(length(G))

                print("evec[1:100]")
                print(evec[1:100])
                # evec = t(I_mat) %*% evec
                # XVsample0_e = XVsample0_e_list[[ne]]
                GE <- G0 * evec
                # GE = G * evec
                # GE = G0
                print("length(GE)")
                print(length(GE))
                GE_tilde <- GE - obj.noK$XXVX_inv %*% (obj.noK$XV %*% GE)
                # GE_tilde_new = GE - I_mat %*% XXVXsample_inv0 %*%  (XVsample0_e %*% G0sample)
                # print("sum(GE_tilde != GE_tilde_new)")
                # print(sum((GE_tilde -  GE_tilde_new)^2))

                # print(obj.noK$XXVX_inv[1:2,])
                # XXVX_invtemp = I_mat %*% XXVXsample_inv0
                # print(XXVX_invtemp[1:2,])

                # obj.glmm.null$eMat[,ne] = GE_tilde
                getildeMat <- cbind(getildeMat, GE_tilde)

                Sigma_iGE <- getSigma_G_multiV(W, tauVecNew, GE_tilde, maxiterPCG, tolPCG, LOCO = FALSE)
                var1GE <- crossprod(GE_tilde, Sigma_iGE) - crossprod(GE_tilde, Sigma_iX) %*% solve(crossprod(X, Sigma_iX)) %*% crossprod(X, Sigma_iGE)
                var1GE_vec <- c(var1GE_vec, var1GE)
                S_GE <- innerProduct(GE_tilde, obj.glmm.null$residuals * var_weights)
                p_exact_GE <- pchisq(S_GE^2 / var1GE, df = 1, lower.tail = F)
                cat("p_exact_GE ", p_exact_GE, "\n")
                cat("S_GE ", S_GE, "\n")
                cat("var1GE ", var1GE, "\n")

                # I_mat_e = I_mat * evec
                # GE_sample = as.vector(t(G0) %*% I_mat_e)
                # GE_sample_tilde = GE_sample  -  XXVXsample_inv0 %*%  (XVsample0 %*% GE_sample)
                # getilde_sample0_Mat = cbind(getilde_sample0_Mat, GE_sample_tilde)
                if (useSparseGRMforVarRatio) {
                  set_isSparseGRM(useSparseGRMforVarRatio)
                  Sigma_iGE_sparse <- getSigma_G_noV(W, tauVecNew, GE_tilde, maxiterPCG, tolPCG, LOCO = FALSE)
                  var2_a_GE <- crossprod(GE_tilde, Sigma_iGE_sparse)
                  var2sparseGRM_GE <- var2_a_GE[1, 1]
                  var2sparseGE_vec <- c(var2sparseGE_vec, var2sparseGRM_GE)
                } else {
                  if (any(duplicated(obj.glmm.null$sampleID))) {
                    if (useGRMtoFitNULL) {
                      tauVal <- tauVecNew[3]
                    } else {
                      tauVal <- tauVecNew[2]
                    }
                    Sigma_iGE_sparse <- getSigma_G_V(W, tauVal, tauVecNew[1], GE_tilde, maxiterPCG, tolPCG)
                    var2_a_GE <- crossprod(GE_tilde, Sigma_iGE_sparse)
                    var2sparseGRM_GE <- var2_a_GE[1, 1]
                    var2sparseGE_vec <- c(var2sparseGE_vec, var2sparseGRM_GE)
                  } else {
                    var2sparseGE_vec <- c(var2sparseGE_vec, var1GE)
                  }
                }
              }
            }


            G_noXadj <- as.vector(G0sample - mean(G0sample))
            G0_sample_tilde <- G0sample - XXVXsample_inv0 %*% XVsample0 %*% G0sample

            # if(useSparseGRMforVarRatio){
            # 	set_isSparseGRM(useSparseGRMforVarRatio)
            # 	Sigma_iG = getSigma_G_noV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO=FALSE)
            # 	var2_a = t(G) %*% Sigma_iG
            # 	var2sparseGRM = var2_a[1,1]
            # 	cat("var2sparseGRM ", var2sparseGRM, "\n")
            # 	varRatio_sparseGRM_vec = c(varRatio_sparseGRM_vec, var1/var2sparseGRM)

            # }else{
            if (any(duplicated(obj.glmm.null$sampleID))) {
              if (useGRMtoFitNULL) {
                tauVal <- tauVecNew[3]
              } else {
                tauVal <- tauVecNew[2]
              }
              # Sigma_iG = getSigma_G_V(W, tauVal, tauVecNew[1], G, maxiterPCG, tolPCG)
              # Sigma_iG = getSigma_G_multiV(W, tauVal, tauVecNew[1], G, maxiterPCG, tolPCG)
              # var2_a = t(G) %*% Sigma_iG
              if (isStoreSigma) {
                Sigma_iG <- (obj.glmm.null$spSigma) %*% G0_sample_tilde
                var2_a <- crossprod(G0_sample_tilde, Sigma_iG)
              } else {
                G0_sample_tilde_I <- as.vector(I_mat %*% G0_sample_tilde)
                Sigma_iG <- getSigma_G_multiV(W, tauVecNew, G0_sample_tilde_I, maxiterPCG, tolPCG, LOCO = FALSE)
                var2_a <- crossprod(G0_sample_tilde_I, Sigma_iG)
              }
              var2sparseGRM <- var2_a[1, 1]
              cat("var2sparseGRM Here ", var2sparseGRM, "\n")
              varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, var1 / var2sparseGRM)
            } else {
              varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, 1)
            }
            # }


            if (obj.glmm.null$traitType == "binary") {
              tmpW <- mu * (1 - mu) * var_weights
              tmpWI <- as.vector(crossprod(tmpW, I_mat))
              var2null <- innerProduct(tmpW, G * G)
              var2null_sample <- innerProduct(tmpWI, G0_sample_tilde * G0_sample_tilde)
              var2null_noXadj <- innerProduct(tmpWI, G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  var22nullGE <- innerProduct(tmpW, GE_tilde * GE_tilde)
                  # var22nullGE = innerProduct(as.vector(t(mu*(1-mu)*var_weights) %*% I_mat), as.vector(GE_sample_tilde*GE_sample_tilde))


                  # var22nullGE = innerProduct(mu*(1-mu)*var_weights, GE_tilde*GE_tilde)
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "quantitative") {
              var2null <- innerProduct(G, G * var_weights)
              var2null_sample <- innerProduct(G0_sample_tilde, G0_sample_tilde * var_weights)
              var2null_noXadj <- innerProduct(G_noXadj, G_noXadj * as.vector(crossprod(var_weights, I_mat)))
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  var22nullGE <- innerProduct(GE_tilde, GE_tilde * var_weights)
                  # var22nullGE = innerProduct(as.vector(GE_sample_tilde), as.vector(GE_sample_tilde)*as.vector(t(var_weights) %*% I_mat))
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "count") {
              tmpW <- mu * var_weights
              tmpWI <- as.vector(crossprod(tmpW, I_mat))
              var2null <- innerProduct(tmpW, G * G)
              # cat("mean(G0_sample_tilde) ", mean(G0_sample_tilde), "\n")
              # var2null_new = innerProduct(as.vector(t(mu*var_weights) %*% I_mat), G0_sample_tilde*G0_sample_tilde)
              var2null_sample <- innerProduct(tmpWI, G0_sample_tilde * G0_sample_tilde)
              # cat("var2null ", var2null, "\n")
              # cat("var2null_new ", var2null_new, "\n")
              # muI <- as.vector(crossprod(mu, I_mat)) * as.vector(var_weights)
              var2null_noXadj <- innerProduct(tmpWI, G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  GE_tilde <- getildeMat[, ne]
                  # GE_sample_tilde = getilde_sample0_Mat[,ne]
                  # cat("mean(GE_sample_tilde) ", mean(GE_sample_tilde), "\n")
                  var22nullGE <- innerProduct(tmpW, GE_tilde * GE_tilde)
                  # var22nullGE = innerProduct(as.vector(t(mu*var_weights) %*% I_mat), as.vector(GE_sample_tilde*GE_sample_tilde))
                  # cat("var22nullGE_old ", var22nullGE_old, "\n")
                  # cat("var22nullGE ", var22nullGE, "\n")

                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            } else if (obj.glmm.null$traitType == "count_nb") {
              var2null <- innerProduct(W, G * G) ## To update
              var2null_sample <- innerProduct(as.vector(crossprod(W, I_mat)), G0_sample_tilde * G0_sample_tilde)
              var2null_noXadj <- innerProduct(as.vector(crossprod(W, I_mat)), G_noXadj * G_noXadj)
              var2nullGE_vec <- NULL
              if (!is.null(obj.glmm.null$eMat)) {
                for (ne in 1:ncol(obj.glmm.null$eMat)) {
                  # GE_tilde = getildeMat[,ne]
                  GE_sample_tilde <- getilde_sample0_Mat[, ne]
                  # var22nullGE = innerProduct(W, GE_tilde*GE_tilde)
                  var22nullGE <- innerProduct(as.vector(crossprod(W, I_mat)), as.vector(GE_sample_tilde * GE_sample_tilde))
                  var2nullGE_vec <- c(var2nullGE_vec, var22nullGE)
                }
              }
            }

            cat("mu\n")
            print(mu[1:100])
            cat("AC ", AC, "\n")
            # cat("S ", S*sqrt(AC), "\n")
            cat("var1 ", var1, "\n")
            cat("var2null ", var2null, "\n")
            cat("var2null_noXadj ", var2null_noXadj, "\n")
            # cat("p_approx ", p_approx, "\n")
            # cat("p_approx_true ", p_approx_true, "\n")
            varRatio_NULL_vec <- c(varRatio_NULL_vec, var1 / var2null)
            varRatio_NULL_sample_vec <- c(varRatio_NULL_sample_vec, var1 / var2null_sample)
            varRatio_NULL_noXadj_vec <- c(varRatio_NULL_noXadj_vec, var1 / var2null_noXadj)
            if (!is.null(obj.glmm.null$eMat)) {
              varRatio_NULL_eg_mat <- rbind(varRatio_NULL_eg_mat, var1GE_vec / var2nullGE_vec)
              varRatio_sparse_eg_mat <- rbind(varRatio_sparse_eg_mat, var1GE_vec / var2sparseGE_vec)
            }
            # indexInMarkerList = indexInMarkerList + 1
            numTestedMarker <- numTestedMarker + 1
          } else {
            indexInMarkerList <- indexInMarkerList + 1
          }

          if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
            numTestedMarker <- numMarkers0
          }
        } # end of while(numTestedMarker < numMarkers)


        print("varRatio_NULL_vec")
        print(varRatio_NULL_vec)
        print("varRatio_NULL_noXadj_vec")
        print(varRatio_NULL_noXadj_vec)

        # ratioCV = calCV(varRatio_NULL_vec)
        ratioCV <- calCV(varRatio_NULL_noXadj_vec)

        if (ratioCV > ratioCVcutoff) {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
          numMarkers0 <- numMarkers0 + 10
          cat("try ", numMarkers0, " markers\n")
        } else {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
        }

        if (indexInMarkerList - 1 == length(listOfMarkersForVarRatio[[k]])) {
          ratioCV <- ratioCVcutoff
          cat("no more markers are available in the MAC category ", k, "\n")
          print(indexInMarkerList - 1)
        }
      } # end of while(ratioCV > ratioCVcutoff)

      if (length(varRatio_sparseGRM_vec) > 0) {
        cat("varRatio_sparseGRM_vec\n")
        print(varRatio_sparseGRM_vec)

        varRatio_sparse <- mean(varRatio_sparseGRM_vec)
        cat("varRatio_sparse", varRatio_sparse, "\n")
        varRatioTable <- rbind(varRatioTable, c(varRatio_sparse, "sparse", k))
      }
      varRatio_null <- mean(varRatio_NULL_vec)
      varRatio_null_sample <- mean(varRatio_NULL_sample_vec)
      cat("varRatio_null", varRatio_null, "\n")

      varRatio_null_noXadj <- mean(varRatio_NULL_noXadj_vec)
      cat("varRatio_null_noXadj", varRatio_null_noXadj, "\n")

      if (!is.null(obj.glmm.null$eMat)) {
        varRatio_NULL_eg_vec <- as.vector(colMeans(varRatio_NULL_eg_mat))
        varRatio_sparse_eg_vec <- as.vector(colMeans(varRatio_sparse_eg_mat))
        for (ne in 1:ncol(obj.glmm.null$eMat)) {
          varRatioTable <- rbind(varRatioTable, c(varRatio_NULL_eg_vec[ne], "null", 0))
          varRatioTable <- rbind(varRatioTable, c(varRatio_sparse_eg_vec[ne], "sparse", 0))
        }
        print(varRatio_NULL_eg_vec)
        print(varRatio_NULL_eg_mat)
      }

      varRatioTable <- rbind(varRatioTable, c(varRatio_null, "null", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_noXadj, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(varRatio_null_sample, "null_sample", k))


      # varRatioTable = rbind(varRatioTable, c(varRatio_null_noXadj, "null", k))
    } else { # if(cateVarRatioVec[k] == 1)
      varRatioTable <- rbind(varRatioTable, c(1, "null", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_noXadj", k))
      varRatioTable <- rbind(varRatioTable, c(1, "null_sample", k))
      if (length(varRatio_sparseGRM_vec) > 0) {
        varRatioTable <- rbind(varRatioTable, c(1, "sparse", k))
      }
    }
  } # for(k in 1:length(listOfMarkersForVarRatio)){
  write.table(varRatioTable, varRatioOutFile, quote = F, col.names = F, row.names = F)
  data <- read.table(varRatioOutFile, header = F)
  print(data)
}


# Fits the null glmm
glmmkin.ai_PCG_Rcpp_multiV <- function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau = c(0, 0), fixtau = c(0, 0), maxiter = 20, tol = 0.02, verbose = TRUE, nrun = 30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, covarianceIdxMat = NULL, isStoreSigma = FALSE, useSparseGRMtoFitNULL = TRUE, useGRMtoFitNULL = TRUE, isSparseGRMIdentity = FALSE) {
  # Fits the null generalized linear mixed model for a poisson, binomial, and gaussian
  # Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  # Returns:
  #  model output for the null glmm

  t_begin <- proc.time()
  print(t_begin)
  subSampleInGeno <- subPheno$IndexGeno
  if (is.null(subPheno$IndexGeno)) {
    subSampleInGeno <- subPheno$IndexPheno
  }
  if (verbose) {
    print("Start reading genotype plink file here")
  }

  # print("subSampleInGeno")
  # print(subSampleInGeno)

  print(head(subPheno))
  set_dup_sample_index_inR(as.numeric(factor(subPheno$IID, levels = unique(subPheno$IID))))


  print("length(indicatorGenoSamplesWithPheno)")
  print(length(indicatorGenoSamplesWithPheno))
  # if((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | (skipVarianceRatioEstimation)){
  if (bedFile != "" & !useSparseGRMtoFitNULL & useGRMtoFitNULL) {
    print("HEREHRE")

    re1 <- system.time({
      setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
    })
  }
  if (verbose) {
    print("Genotype reading is done")
  }

  if (LOCO) {
    MsubIndVec <- getQCdMarkerIndex()
    chrVec <- data.table::fread(bimFile, header = F)[, 1]
    chrVec <- chrVec[which(MsubIndVec == TRUE)]
    updatechrList <- updateChrStartEndIndexVec(chrVec)
    LOCO <- updatechrList$LOCO
    chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
  }

  y <- fit0$y
  n <- length(y)
  X <- model.matrix(fit0)
  offset <- fit0$offset
  if (is.null(offset)) {
    offset <- rep(0, n)
  }

  var_weights <- weights(fit0)


  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu) / mu.eta

  if (is.null(var_weights)) {
    var_weights <- rep(1, length(mu.eta))
  }
  # sqrtW = mu.eta/sqrt(family$variance(mu))
  # }else{
  sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
  # }
  W <- sqrtW^2

  alpha0 <- fit0$coef
  eta0 <- eta


  cat("tauInit")
  print(tauInit)

  tau[1:length(tau)] <- 0
  if (family$family %in% c("poisson", "binomial")) {
    # if(family$family %in% c("binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
    tauInit[1] <- 1
    idxtau <- which(fixtau == 0)
    cat("fixtau ", fixtau, "\n")
    cat("tauInit ", tauInit, "\n")
    if (sum(tauInit[idxtau]) == 0) {
      tau[idxtau] <- 0.1
    } else {
      tau[idxtau] <- tauInit[idxtau]
    }
  } else { #  if(family$family %in% c("poisson", "binomial")) {
    idxtau <- which(fixtau == 0)
    if (sum(tauInit[idxtau]) == 0) {
      tau[1] <- 1
      # tauInit[1] = 1
      tau[idxtau] <- var(Y) / (length(tau))
      # tau[2] = 0
      # tau[2:length(tau)] = 0
      if (abs(var(Y)) < 0.1) {
        stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
      }
    } else {
      tau[fixtau == 0] <- tauInit[fixtau == 0]
    }
  }

  cat("inital tau is ", tau, "\n")

  if (!is.null(covarianceIdxMat)) {
    idxtau2 <- intersect(covarianceIdxMat[, 1], idxtau)
    print("covarianceIdxMat")
    print(covarianceIdxMat)
    print("idxtau2")
    print(idxtau2)
    if (length(idxtau2) > 0) {
      tau[idxtau2] <- 0
    }
    # i_kmat = get_numofV()
    # if(i_kmat > 0){
    Kmatdiag <- getMeanDiagofKmat(LOCO)
    print(Kmatdiag)
    # }
    tau[2:length(tau)] <- tau[2:length(tau)] / Kmatdiag
  }

  print("tau")
  print(tau)


  #### set up weights for variance
  # if(!is.null(var_weights)){
  # 	set_var_weights(var_weights)
  # }


  # if(isStoreSigma){
  #  gen_sp_Sigma_multiV(W, tau)
  # }

  if (isSparseGRMIdentity) {
    tau[2] <- 0
  }

  # print("Here isStoreSigma")
  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  # if(isStoreSigma){
  #  gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  # print("Here isStoreSigma 2")

  re <- getAIScore_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff, LOCO = FALSE)
  tau0 <- tau
  tau0_q2 <- tau[fixtau == 0]
  # tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)
  print("tau0_q2 a")
  print(tau0_q2)
  print("idxtau")
  print(idxtau)


  tau_q2 <- pmax(0, tau0_q2 + tau0_q2^2 * (re$YPAPY - re$Trace) / n)
  tau[idxtau] <- tau_q2

  if (!is.null(covarianceIdxMat)) {
    tau[idxtau[which(idxtau %in% idxtau2)]] <- 0
  }
  print("re$YPAPY")
  print(re$YPAPY)
  print("re$Trace")
  print(re$Trace)

  if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  maxiter_in <- maxiter
  if (isSparseGRMIdentity) {
    tau[2] <- 0
    maxiter_in <- 0
    alpha <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
  }

  for (i in seq_len(maxiter_in)) {
    # W = sqrtW^2

    if (verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
    # use Get_Coef before getAIScore
    t_begin_Get_Coef <- proc.time()
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(W, tau)
    # }
    # cat("eta0 ", eta0, "\n")
    re.coef <- Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)
    t_end_Get_Coef <- proc.time()
    cat("t_end_Get_Coef - t_begin_Get_Coef\n")
    print(t_end_Get_Coef - t_begin_Get_Coef)
    # if(isStoreSigma){
    #  gen_sp_Sigma_multiV(re.coef$W, tau)
    # }


    fit <- fitglmmaiRPCG_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff, LOCO = FALSE)

    t_end_fitglmmaiRPCG <- proc.time()
    cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

    tau <- as.numeric(fit$tau)
    cov <- re.coef$cov
    alpha <- re.coef$alpha
    eta <- re.coef$eta
    Y <- re.coef$Y
    mu <- re.coef$mu

    mu.eta <- family$mu.eta(eta)

    if (is.null(var_weights)) {
      sqrtW <- mu.eta / sqrt(family$variance(mu))
    } else {
      sqrtW <- mu.eta / sqrt(1 / as.vector(var_weights) * family$variance(mu))
    }
    W <- sqrtW^2


    print(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol))
    cat("tau: ", tau, "\n")
    cat("tau0: ", tau0, "\n")

    # if(family$family == "gaussian"){
    # if(tau[1]<=0){
    #  tau[1] = tau[1] + 0.1
    #  #stop("ERROR! The first variance component parameter estimate is 0\n")
    # }
    # }

    # if(sum(tau[2:length(tau)]) == 0) break
    # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    if (sum(tau[2:length(tau)]) == 0) {
      break
      # tau[2:length(tau)] = rep(0.1,length(tau)-1)
    } else {
      if (max(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)) < tol) break

      if (max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i <- maxiter
        break
      }
    }
  }

  if (verbose) cat("\nFinal ", tau, ":\n")


  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(W, tau)
  # }

  # added these steps after tau is estimated 04-14-2018

  re.coef <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = FALSE, var_weights = var_weights)

  cov <- re.coef$cov
  alpha <- re.coef$alpha
  eta <- re.coef$eta
  Y <- re.coef$Y
  mu <- re.coef$mu
  converged <- ifelse(i < maxiter, TRUE, FALSE)

  # var_weights = NULL
  # if(!is.null(var_weights)){
  #  res = (y - mu) * sqrt(var_weights)
  # }else{
  res <- y - mu
  # }

  if (family$family == "binomial") {
    mu2 <- mu * (1 - mu)
    traitType <- "binary"
  } else if (family$family == "poisson") {
    mu2 <- mu
    traitType <- "count"
  } else if (family$family == "gaussian") {
    mu2 <- rep((1 / (tau[1])), length(res))
    # mu2 = rep(1,length(res))
    traitType <- "quantitative"
  }

  # if(isCovariateTransform & hasCovariate){
  # if(!is.null(out.transform) & is.null(fit0$offset)){
  #  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  # }else{
  #  coef.alpha = alpha
  # }

  # mu2 = mu * (1-mu)

  # if(!is.null(var_weights)){
  mu2_rescaled <- mu2 * var_weights
  y_rescaled <- y * var_weights
  mu_rescaled <- mu * var_weights
  # }else{


  # }

  if (!isCovariateOffset) {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, barcode = subPheno$barcode, obj.noK = obj.noK, y = y, X = X, traitType = traitType, isCovariateOffset = isCovariateOffset)
  } else {
    obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, barcode = subPheno$barcode, obj.noK = obj.noK, y = y, X = Xorig, traitType = traitType, isCovariateOffset = isCovariateOffset)
  }

  glmmResult$varWeights <- var_weights
  # LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme
  # glmmResult$Sigma_iX = re.coef$Sigma_iX

  glmmResult$LOCO <- LOCO
  t_end_null <- proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  # if(isStoreSigma){
  #   gen_sp_Sigma_multiV(re.coef$W, tau)
  # }
  if (!isLowMemLOCO & LOCO) {
    # if(isStoreSigma){
    #   gen_sp_Sigma_multiV(re.coef$W, tau)
    # }
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult <- list()
    for (j in 1:22) {
      startIndex <- chromosomeStartIndexVec[j]
      endIndex <- chromosomeEndIndexVec[j]
      if (!is.na(startIndex) && !is.na(endIndex)) {
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j - 1)
        t_begin_Get_Coef_LOCO <- proc.time()
        re.coef_LOCO <- Get_Coef_multiV(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, LOCO = TRUE, var_weights = varWeights)
        t_end_Get_Coef_LOCO <- proc.time()
        cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov <- re.coef_LOCO$cov
        alpha <- re.coef_LOCO$alpha
        eta <- re.coef_LOCO$eta
        Y <- re.coef_LOCO$Y
        mu <- re.coef_LOCO$mu
        # mu2 = mu * (1-mu)
        # mu2 = mu

        res <- y - mu


        if (family$family == "binomial") {
          mu2 <- mu * (1 - mu)
        } else if (family$family == "poisson") {
          mu2 <- mu
        } else if (family$family == "gaussian") {
          mu2 <- rep((1 / (tau[1])), length(res))
        }


        if (!is.null(out.transform) & is.null(fit0$offset)) {
          coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        } else {
          coef.alpha <- alpha
        }

        mu2_rescaled <- mu2 * var_weights
        mu_rescaled <- mu * var_weights


        if (!isCovariateOffset) {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X)
        } else {
          obj.noK <- ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig)
        }
        glmmResult$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)
      } else {
        glmmResult$LOCOResult[[j]] <- list(isLOCO = FALSE)
      }
    }
  }

  if (isLowMemLOCO & LOCO) {
    glmmResult$chromosomeStartIndexVec <- chromosomeStartIndexVec
    glmmResult$chromosomeEndIndexVec <- chromosomeEndIndexVec
  }
  return(glmmResult)
}
