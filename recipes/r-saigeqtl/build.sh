#!/bin/bash

$R CMD INSTALL --build .

pushd extdata
install "step1_fitNULLGLMM_qtl.R" "step2_tests_qtl.R" "step3_gene_pvalue_qtl.R" "makeGroupFile.R" "createSparseGRM.R" "${PREFIX}/bin"
popd
