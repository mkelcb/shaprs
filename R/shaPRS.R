#' Create lFDR corrected Q-test statistics for each SNP
#'
#' it performs:
#' (1) modified Cochran's Q-test which optionally adjusts for overlapping controls
#' (2) lFDR estimation on the p-values from the above Cochran's Q-test
#' (3) lists SNPs that are fail the heterogeneity test at specified thresholds (optional)
#'
#' @param inputData summary statistics table that has header with the following columns: SNP	CHR	BP	Beta_A	SE_A	Beta_B	SE_B
#' @param rho estimate of correlation between studies due to shared subjects. 0 for no overlap and 1 for complete overlap. default: 0. Obtain this from shaPRS_rho()
#' @param thresholds vector of thresholds to be used to create list of SNPs (default empty)
#' @return returns object with two fields, (1) lFDRTable: a 2 column file with the following signature SNPID lFDR (2) hardThresholds list of SNPids that failed the heterogeneity test at each threshold
#'
#' @importFrom stats na.omit pchisq pnorm cor
#'
#' @examples
#' inputDataLoc <- system.file("extdata", "shapersToydata.txt", package = "shaPRS")
#' inputData= read.table(inputDataLoc, header = TRUE)
#' results = shaPRS_adjust(inputData, thresholds=c(0.5,0.99))
#'
#' @export
shaPRS_adjust = function(inputData, rho = 0, thresholds =  vector()) {

  # 1. Cochran's Q-test formula: from 'Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects' ncbi.nlm.nih.gov/pmc/articles/PMC2790578/
  Vhat = (inputData$SE_A^2 + inputData$SE_B^2 - 2 * rho * inputData$SE_A * inputData$SE_B)
  Q_adjusted = (inputData$Beta_A - inputData$Beta_B)^2/ Vhat
  df=2-1 # degrees of freedom, 2 studies -1
  Q_adjusted_pvals = pchisq(Q_adjusted, df = df, lower.tail = F)

  # 2. lFDR estimation
  lfdr_obj = qvalue::qvalue(p = Q_adjusted_pvals)
  lfdr_qvals <- lfdr_obj$lfdr

  # 3. prepare table of lFDR values for each SNP
  lFDRTable <- data.frame(inputData$SNP, lfdr_qvals)

  # 4. Create list for each threshold of SNPs that fail the heterogeneity test at the specified thresholds
  hard_threshold_results = list()
  if (length(thresholds) > 0) {
    for (i in 1:length(thresholds)) {
      currentThreshold=thresholds[i]
      sigSNPs = lFDRTable[which(lFDRTable$lfdr_qvals <currentThreshold),]

      # SNPs which are heterogeneous, will be sourced from the subphenotpye
      CompositePheno = inputData # start from the composite pheno
      hard_threshold_results[[i]] = sigSNPs$inputData.SNP
    }
  }



  results <- list("lFDRTable" = lFDRTable, "hardThresholds" = hard_threshold_results)
  return(results)
}


#' Composite shaPRS: produce summary statistics according to hard thresholds
#'
#' Loops through a list of hard thresholds, and swaps between the Combined and Subphenotype estimates
#' based on the list of SNPs found to be heterogenous for a given threshold. This process would swap in
#' the subphenotype estimates for the SNPs found to be heterogeneous.
#'
#' @param subpheno  Subphenotype LDPred formatted GWAS summary statistics table  that has header with the following columns: chr	pos	SNP	A1	A2	Freq1.Hapmap	b	se	p	N
#' @param CombinedPheno Combined phenotype LDPred formatted GWAS summary statistics table, same signature
#' @param hard_threshold_results a list of hard thresholds, in the format produced by shaPRS_adjust()
#' @return returns a list of LDPred formatted summary statistics tables, one for each threshold
#'
#' @examples
#' inputDataLoc <- system.file("extdata", "shapersToydata.txt", package = "shaPRS")
#' inputData= read.table(inputDataLoc, header = TRUE)
#' results = shaPRS_adjust(inputData, thresholds=c(0.5,0.99))
#' subphenoLoc <- system.file("extdata", "phenoA_sumstats", package = "shaPRS")
#' CombinedPhenoLoc <- system.file("extdata", "Combined_sumstats", package = "shaPRS")
#' subpheno= read.table(subphenoLoc, header = TRUE)
#' CombinedPheno= read.table(CombinedPhenoLoc, header = TRUE)
#' CompositeSumstats = shaPRS_composite(subpheno, CombinedPheno, results$hardThresholds)
#'
#' @export
shaPRS_composite = function(subpheno, CombinedPheno, hard_threshold_results) {

  allCompositeResults= list()
  for (i in 1:length(hard_threshold_results)) {
    heterogeneous_SNPs = hard_threshold_results[[i]]

    # SNPs which are heterogeneous, will be sourced from the subphenotpye
    CompositePheno = CombinedPheno # start by duplicating the combined

    subphenoIndices = match(heterogeneous_SNPs,subpheno$SNP) # find the indices of the SNPs to be swapped
    swappedSNPs = subpheno[subphenoIndices,]
    CompositePheno[subphenoIndices,] = subpheno[subphenoIndices,]

    allCompositeResults[[i]] = CompositePheno
  }

  return(allCompositeResults)
}



#' Blended shaPRS: produce summary statistics according to a continuous weighting scheme
#'
#' This function continuously blends the sub-phenotype and the combined phenotype summary statistics
#' and generates an LDPred formatted table.
#'
#' @param subpheno  Subphenotype LDPred formatted GWAS summary statistics table  that has header with the following columns: chr	pos	SNP	A1	A2	Freq1.Hapmap	b	se	p	N
#' @param CombinedPheno Combined phenotype LDPred formatted GWAS summary statistics table, same signature
#' @param blendingFactors a headerless 2 column table of the SNPid and lFDR, (produced by shaPRS_adjust)
#' @return returns an LDPred formatted summary statistics table
#'
#' @importFrom stats na.omit pchisq pnorm cor
#'
#' @examples
#' subphenoLoc <- system.file("extdata", "phenoA_sumstats", package = "shaPRS")
#' CombinedPhenoLoc <- system.file("extdata", "Combined_sumstats", package = "shaPRS")
#' blendFactorLoc <- system.file("extdata", "myOutput_SNP_lFDR", package = "shaPRS")
#' subpheno= read.table(subphenoLoc, header = TRUE)
#' CombinedPheno= read.table(CombinedPhenoLoc, header = TRUE)
#' blendingFactors= read.table(blendFactorLoc, header = FALSE)
#' blendedSumstats = shaPRS_blend(subpheno, CombinedPheno, blendingFactors)
#'
#' @export
shaPRS_blend = function(subpheno, CombinedPheno, blendingFactors) {
  # 1. merge the 3 tables together by RSid, so they are always aligned, x = subpheno  and    y = CombinedPheno
  subpheno_CombinedPheno = merge(subpheno,CombinedPheno,by.x = "SNP",by.y = "SNP")
  subpheno_CombinedPheno_blending = merge(subpheno_CombinedPheno,blendingFactors,by.x = "SNP",by.y = "V1")

  # 2. Blend the subpheno and CombinedPheno together and create new summary statistics via following logic:
  # blendingfactor = lFDR = prob that it is the null, that it was CombinedPheno NOT the subpheno
  # basic formula: beta_phenoA x P(phenoA != phenoB) + beta_CombinedPheno x P(phenoA == phenoB)
  # Var(aX+bY)=a^2 var(X) + b^2 var(Y)  +2abcov(X,Y)
  # se is sqrt var, so:
  # sqrt( Var(aX+bY) )= sqrt( a^2*var(X) + b^2*var(Y) +2abcov(X,Y) )
  blendedSE = sqrt(  (1-subpheno_CombinedPheno_blending$V2)^2 * subpheno_CombinedPheno_blending$se.x^2 +  subpheno_CombinedPheno_blending$se.y^2 * subpheno_CombinedPheno_blending$V2^2 + 2*subpheno_CombinedPheno_blending$V2*(1-subpheno_CombinedPheno_blending$V2)*cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)* subpheno_CombinedPheno_blending$se.x*subpheno_CombinedPheno_blending$se.y )
  blendedBeta=subpheno_CombinedPheno_blending$b.x * (1-subpheno_CombinedPheno_blending$V2) + subpheno_CombinedPheno_blending$b.y * subpheno_CombinedPheno_blending$V2
  blendedp=2*pnorm( abs(blendedBeta)/blendedSE,lower.tail=FALSE)

  # 3. create new data frame to store the new summary stats
  blendedSumstats = data.frame(subpheno_CombinedPheno_blending$chr.x, subpheno_CombinedPheno_blending$pos.x, subpheno_CombinedPheno_blending$SNP, subpheno_CombinedPheno_blending$A1.x, subpheno_CombinedPheno_blending$A2.x, subpheno_CombinedPheno_blending$Freq1.Hapmap.x,
                               blendedBeta,
                               blendedSE,
                               blendedp,
                               round(subpheno_CombinedPheno_blending$N.x * (1-subpheno_CombinedPheno_blending$V2) + subpheno_CombinedPheno_blending$N.y * subpheno_CombinedPheno_blending$V2)
  )
  colnames(blendedSumstats) = colnames(CombinedPheno)
  blendedSumstats= blendedSumstats[match(subpheno$SNP, blendedSumstats$SNP),]
  blendedSumstats <- na.omit(blendedSumstats) # remove any NAs of SNPs that couldn't be matched

  return(blendedSumstats)
}

#' Calculate rho to be used in shaPRS_adjust()
#'
#' Convenience function to estimate the correlation between  two studies with overlapping controls
#' for more details see:
#' Lin et al. 'Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects' (2009
#' ncbi.nlm.nih.gov/pmc/articles/PMC2790578
#'
#' @param nkl0 number of controls overlapping between studies
#' @param nk1 number of cases in study k
#' @param nk0 number of controls in study k
#' @param nl1 number of cases in study l
#' @param nl0 number of controls in study l
#' @return returns real value of the approximate correlation
#'
#' @examples
#' rho = shaPRS_rho(nkl0 = 9492,nk1 = 3810, nk0= 9492, nl1= 3765,nl0= 9492)
#'
#' @export
shaPRS_rho = function(nkl0,nk1, nk0, nl1,nl0) {
  nk=nk1+nk0 # total number indis in k
  nl= nl1+nl0 # total number of indis in l
  approx_cor= (nkl0 * sqrt(nk1*nl1 / (nk0*nl0) )   ) / sqrt(nk*nl)
  return(approx_cor )
}



