#' Create lFDR corrected Q-test statistics for each SNP
#'
#' it performs:
#' (1) modified Cochran's Q-test which optionally adjusts for overlapping controls
#' (2) lFDR estimation on the p-values from the above Cochran's Q-test
#' (3) lists SNPs that fail the heterogeneity test at specified thresholds (optional)
#'
#' @param inputData summary statistics table that has header with the following columns: SNP	CHR	BP	Beta_A	SE_A	A1.x	A2.x	Beta_B	SE_B	A1.y	A2.y
#' @param rho estimate of correlation between studies due to shared subjects. 0 for no overlap and 1 for complete overlap. default: 0. Obtain this from shaPRS_rho()
#' @param thresholds vector of thresholds to be used to create list of SNPs (default empty)
#' @return returns object with two fields, (1) lFDRTable: a 3 column file with the following signature SNP lFDR Qval (2) hardThresholds list of SNPids that failed the heterogeneity test at each threshold
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

  # remove non numeric data from cols that must be numeric (must cast to 'characer' otherwise R may turn a value eg 0.8249 into a large integer like 8150 on nix systems)
  inputData <- inputData[!is.na(as.numeric(as.character(inputData$SE_A))),]
  inputData <- inputData[!is.na(as.numeric(as.character(inputData$SE_B))),]
  inputData <- inputData[!is.na(as.numeric(as.character(inputData$Beta_A))),]
  inputData <- inputData[!is.na(as.numeric(as.character(inputData$Beta_B))),]
  # now actually cast them to numeric
  inputData$SE_A = as.numeric(as.character(inputData$SE_A ))
  inputData$SE_B = as.numeric(as.character(inputData$SE_B ))
  inputData$Beta_A = as.numeric(as.character(inputData$Beta_A ))
  inputData$Beta_B = as.numeric(as.character(inputData$Beta_B ))

  inputData = alignStrands(inputData)


  # 0. Reverse effect sizes alleles
  misalignedAlleleIndices = which( as.character(inputData$A1.x) != as.character(inputData$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
  inputData$Beta_B[misalignedAlleleIndices] = -inputData$Beta_B[misalignedAlleleIndices] # flip effects
  if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )

  # 1. Cochran's Q-test formula: from 'Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects' ncbi.nlm.nih.gov/pmc/articles/PMC2790578/
  Vhat = (inputData$SE_A^2 + inputData$SE_B^2 - 2 * rho * inputData$SE_A * inputData$SE_B)
  Q_vals = (inputData$Beta_A - inputData$Beta_B)^2/ Vhat
  df=2-1 # degrees of freedom, 2 studies -1
  Q_pvals = pchisq(Q_vals, df = df, lower.tail = F)

  # 2. lFDR estimation
  lfdr_obj = qvalue::qvalue_truncp(p = Q_pvals)
  lfdr_qvals <- lfdr_obj$lfdr

  # 3. prepare table of lFDR values for each SNP
  lFDRTable <- data.frame(inputData$SNP, lfdr_qvals, Q_vals)
  colnames(lFDRTable) = c("SNP", "lFDR", "Qval")
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


#alleles = c("G", "G","C","G")
flipStrand = function(alleles) {
  allelesFlipped = alleles
  whereAsare = which(alleles == "A")
  whereGsare = which(alleles == "G")
  allelesFlipped[which(allelesFlipped == "T")] = "A" # flip Ts to As
  allelesFlipped[which(allelesFlipped == "C")] = "G" # flip Cs to Gs
  allelesFlipped[whereAsare] = "T" # flip As to Ts
  allelesFlipped[whereGsare] = "C" # flip Gs to Cs
  return(allelesFlipped)
}


# inputData = subphenoA_B
#inputData = subpheno_otherPheno_blending
#inputData = sumstatsDataAll

#' Aligns strands between 2 summary statistics data
#'
#'
#'
#' @param inputData dataframe of both studies with columns:  A1.x, A2.x, A1.y, A2.y
#' @param A1.x column for A1 allele in study 1
#' @param A2.x column for A2 allele in study 1
#' @param A1.y column for A1 allele in study 2
#' @param A2.y column for A1 allele in study 2
#' @return returns real value of the approximate correlation
#'
#'
#' @export
alignStrands = function(inputData, A1.x ="A1.x", A2.x ="A2.x", A1.y ="A1.y", A2.y ="A2.y") {

  # exclude ambiguous SNPs
  ambiguousSNPIndices = which(inputData[,A1.x] == "G" & inputData[,A2.x] == "C" |  inputData[,A1.x] == "C" & inputData[,A2.x] == "G" | inputData[,A1.x] == "A" & inputData[,A2.x] == "T"  | inputData[,A1.x == "T"] &  inputData[,A2.x] ==  "A" |
                              inputData[,A1.y] == "G" & inputData[,A2.y ]== "C" |  inputData[,A1.y] == "C" & inputData[,A2.y] == "G" | inputData[,A1.y] == "A" & inputData[,A2.y] == "T"  | inputData[,A1.y == "T"] &  inputData[,A2.y] ==  "A")

  print(paste0("removed ", length(ambiguousSNPIndices), " ambiguous SNPs out of ", nrow(inputData), " variants"))
  if(length(ambiguousSNPIndices) > 0) inputData = inputData[-ambiguousSNPIndices,] # otherwise R would remove all as R is shit

  # 0. Reverse strands: (this assumes that ambigous SNPs have been excluded prior to this)
  # cache flipped strands for phenoB
  inputData$A1.y_flipped =  flipStrand(inputData[,A1.y])
  inputData$A2.y_flipped =  flipStrand(inputData[,A2.y])
  # if 2nd sumstats CAN be matched to 1st, IE
  #                                                                                 regular match                                                                                                                         reverse match                                                                                                    flipped match                                                                                                    reverse flipped match
  matchedIndices = which(as.character(inputData[,A1.x]) == as.character(inputData[,A1.y]) & as.character(inputData[,A2.x]) == as.character(inputData[,A2.y]) |  as.character(inputData[,A1.x]) == as.character(inputData[,A2.y]) & as.character(inputData[,A2.x]) == as.character(inputData[,A1.y]) | as.character(inputData[,A1.x]) == as.character(inputData$A1.y_flipped) & as.character(inputData[,A2.x]) == as.character(inputData$A2.y_flipped)  | as.character( inputData[,A1.x]) == as.character(inputData$A2.y_flipped) &  as.character(inputData[,A2.x]) == as.character(inputData$A1.y_flipped) )
  print(paste0("matched ", length(matchedIndices), " out of ", nrow(inputData), " variants"))
  # as there may be non-SNPs, we need to cast them as character


  # exclude unmatchables
  inputData = inputData[matchedIndices,]
  #inputData = inputData_orig


  # flip those A1/A2s which were flipped matches
  #                                                                                       flipped match                                                                                                                                reverse flipped match
  flippedIndices = which(as.character(inputData[,A1.x]) == as.character(inputData$A1.y_flipped) & as.character(inputData[,A2.x]) == as.character(inputData$A2.y_flipped) | as.character(inputData[,A1.x]) == as.character(inputData$A2.y_flipped) & as.character(inputData[,A2.x]) ==  as.character(inputData$A1.y_flipped) )
 # print( head(inputData[flippedIndices,]) )

  print(paste0("flipped strand for ", length(flippedIndices), " variants"))
  inputData[flippedIndices,A1.y] = as.character(inputData$A1.y_flipped[flippedIndices])
  inputData[flippedIndices,A2.y] = as.character(inputData$A2.y_flipped[flippedIndices])

  return(inputData)
}



#' Blended shaPRS (with overlapping datasets): produce summary statistics according to a continuous weighting scheme
#'
#' This function continuously blends the two sub-phenotype statistics
#' and generates an LDPred formatted table.
#'
#' @param subpheno  Subphenotype LDPred formatted GWAS summary statistics table  that has header with the following columns: chr	pos	SNP	A1	A2	Freq1.Hapmap	b	se	p	N
#' @param subpheno_other dataframe for other subpheno of the same signature
#' @param blendingFactors a 3 column table of: SNP lFDR Qval, (produced by shaPRS_adjust)
#' @param rho (optional) overlap between studies
#' @return returns an LDPred formatted summary statistics table
#'
#' @importFrom stats na.omit pchisq pnorm cor
#'
#' @examples
#' subphenoLoc <- system.file("extdata", "phenoA_sumstats", package = "shaPRS")
#' subpheno_otherLoc <- system.file("extdata", "phenoB_sumstats", package = "shaPRS")
#' blendFactorLoc <- system.file("extdata", "myOutput_SNP_lFDR", package = "shaPRS")
#' subpheno= read.table(subphenoLoc, header = TRUE)
#' subpheno_other= read.table(subpheno_otherLoc, header = TRUE)
#' blendingFactors= read.table(blendFactorLoc, header = TRUE)
#' blendedSumstats = shaPRS_blend_overlap(subpheno, subpheno_other, blendingFactors)
#'
#' @export
shaPRS_blend_overlap = function(subpheno, subpheno_other, blendingFactors, rho = 0) {
  # cast as numeric
  subpheno = RemoveNonNumerics(subpheno)
  subpheno_other = RemoveNonNumerics(subpheno_other)

  # 1.  Merge first the 3 tables together by RSid, so they are always aligned, x = subpheno  and    y = CombinedPheno ( ensure that when we check allele alignment we are comparing the same SNPs
  subpheno_otherPheno = merge(subpheno,subpheno_other,by.x = "SNP",by.y = "SNP")
  subpheno_otherPheno_blending = merge(subpheno_otherPheno,blendingFactors, by.x = "SNP", by.y = "SNP")


  subpheno_otherPheno_blending = alignStrands(subpheno_otherPheno_blending)

  # 2. Align PheB/B alleles
  misalignedAlleleIndices = which( as.character(subpheno_otherPheno_blending$A1.x) != as.character(subpheno_otherPheno_blending$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
  subpheno_otherPheno_blending$b.y[misalignedAlleleIndices] = -subpheno_otherPheno_blending$b.y[misalignedAlleleIndices] # flip effects
  if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )

  # sanitize each input of NAs
  subpheno_otherPheno_blending <- na.omit(subpheno_otherPheno_blending) # remove any NAs of SNPs

  w = subpheno_otherPheno_blending$lFDR
  tao1 = 1/subpheno_otherPheno_blending$se.x^2
  tao2 = 1/subpheno_otherPheno_blending$se.y^2

  # calculate the meta analysis beta coefficients and standard errors
  meta_se = sqrt( (tao1 + tao2 + rho *sqrt(tao1 * tao2) ) / ( (tao1 + tao2)^2 ) ) # when rho ==0, this is identical to sqrt(CovB12), but otherwwise this is different
  meta_coef = (subpheno_otherPheno_blending$b.x*1/subpheno_otherPheno_blending$se.x^2 + subpheno_otherPheno_blending$b.y* 1/subpheno_otherPheno_blending$se.y^2) / (1/subpheno_otherPheno_blending$se.x^2+ 1/subpheno_otherPheno_blending$se.y^2)



  # 3. Blend the subpheno and CombinedPheno together and create new summary statistics via following logic:
  # Theoretical (Chris's orginal):
  CovB12 =  ( 1 + rho * sqrt(tao2/tao1) ) / (tao1 + tao2)
  blendedSE =  sqrt( (1-w)^2 /tao1 + w^2 * (tao1 +tao2 + rho * sqrt(tao1 * tao2) ) / (tao1 + tao2)^2 + 2 * w*(1-w) * CovB12 )

  # Empirical corr
  #CovB12_empirical = cor(subpheno_otherPheno_blending$b.x,meta_coef)* subpheno_otherPheno_blending$se.x*meta_se
  #blendedSE =  sqrt( (1-w)^2 /tao1 + w^2 * (tao1 +tao2 + rho * sqrt(tao1 * tao2) ) / (tao1 + tao2)^2 + 2 * w*(1-w) * CovB12_empirical )



  blendedBeta=subpheno_otherPheno_blending$b.x * (1-subpheno_otherPheno_blending$lFDR) + meta_coef * subpheno_otherPheno_blending$lFDR
  blendedp=2*pnorm( abs(blendedBeta)/blendedSE,lower.tail=FALSE)

  # also need the combined sample size
  CombinedN = subpheno_otherPheno_blending$N.x + subpheno_otherPheno_blending$N.y
  # 3. create new data frame to store the new summary stats
  blendedSumstats = data.frame(subpheno_otherPheno_blending$chr.x, subpheno_otherPheno_blending$pos.x, subpheno_otherPheno_blending$SNP, subpheno_otherPheno_blending$A1.x, subpheno_otherPheno_blending$A2.x, subpheno_otherPheno_blending$Freq1.Hapmap.x,
                               blendedBeta,
                               blendedSE,
                               blendedp,
                               round(subpheno_otherPheno_blending$N.x * (1-subpheno_otherPheno_blending$lFDR) + CombinedN * subpheno_otherPheno_blending$lFDR)
  )
  colnames(blendedSumstats) = colnames(subpheno)
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


#' Generic inverse variance meta analysis
#'
#' Convenience function to produce the combined phenotype estimate if we only have summary stats for phenoA and pheno B
#'
#' @param subpheno dataframe for main subpheno
#' @param subpheno_other dataframe for other subpheno
#' @param rho (optional) overlap between studies
#' @return returns Combinedpheno dataframe that can be plugged into shaPRS_blend or shaPRS_composite
#'
#'
#' @export
inverse_metaAnalaysis = function(subpheno,subpheno_other, rho = 0) {

  # cast as numeric
  subpheno = RemoveNonNumerics(subpheno)
  subpheno_other = RemoveNonNumerics(subpheno_other)

  # 1.  Merge first the tables together by RSid, so they are always aligned, x = subpheno  and    y = subpheno_other ( ensure that when we check allele alignment we are comparing the same SNPs
  subphenoA_B = merge(subpheno,subpheno_other,by.x = "SNP",by.y = "SNP")


  subphenoA_B = alignStrands(subphenoA_B)

  # 2. Align PheB/B alleles
  misalignedAlleleIndices = which( as.character(subphenoA_B$A1.x) != as.character(subphenoA_B$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
  subphenoA_B$b.y[misalignedAlleleIndices] = -subphenoA_B$b.y[misalignedAlleleIndices] # flip effects for phe B
  if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )


  # INVERSE VARIANCE FIXED EFFECT META ANALYSIS: https://en.wikipedia.org/wiki/Inverse-variance_weighting
  meta_coef = (subphenoA_B$b.x*1/subphenoA_B$se.x^2 + subphenoA_B$b.y* 1/subphenoA_B$se.y^2) / (1/subphenoA_B$se.x^2+ 1/subphenoA_B$se.y^2)

  # no overlap, use simple formula
  if(rho == 0) { meta_se = sqrt( (subphenoA_B$se.x^(-2) + subphenoA_B$se.y^(-2) )^(-1) ) }
  else { # there is an overlap, use Chris' updated formula
    tao1 = 1/subphenoA_B$se.x^2
    tao2 = 1/subphenoA_B$se.y^2
    meta_se = sqrt( (tao1 + tao2 + rho *sqrt(tao1 * tao2) ) / ( (tao1 + tao2)^2 ) )
  }


  meta_p = 2*pnorm(-abs(meta_coef/meta_se))


  # 3. create new data frame to store the new summary stats
  CombinedPheno = data.frame(subphenoA_B$chr.x, subphenoA_B$pos.x, subphenoA_B$SNP, subphenoA_B$A1.x, subphenoA_B$A2.x, subphenoA_B$Freq1.Hapmap.x,
                             meta_coef,
                             meta_se,
                             meta_p,
                             subphenoA_B$N.x + subphenoA_B$N.y)

  colnames(CombinedPheno) = colnames(subpheno)
  CombinedPheno= CombinedPheno[match(subpheno$SNP, CombinedPheno$SNP),]
  CombinedPheno <- na.omit(CombinedPheno) # remove any NAs of SNPs that couldn't be matched

  return(CombinedPheno)
}

# helper function that removes and casts all columns that should be numeric as numeric
RemoveNonNumerics = function(subpheno) {
  # remove non numeric data from cols that must be numeric (must cast to 'characer' otherwise R may turn a value eg 0.8249 into a large integer like 8150 on nix systems)
  subpheno <- subpheno[!is.na(as.numeric(as.character(subpheno$b))),]
  subpheno <- subpheno[!is.na(as.numeric(as.character(subpheno$se))),]

  # now actually cast them to numeric
  subpheno$b = as.numeric(as.character(subpheno$b ))
  subpheno$se = as.numeric(as.character(subpheno$se ))
  return(subpheno)
}


#require(compiler)
#enableJIT(3)

#library(Matrix)
#Matrix compiler



#' Generate shaPRS specific LD refernece panel
#'
#' Generates a PRS specific LD reference matrix by blending together two LD ref panels according to
#' shaPRS produced lFDR and standard errors
#'
#' @param pop1LDmatrix LD reference matrix in RDS (dsCMatrix) format for target population
#' @param pop2LDmatrix LD reference matrix in RDS (dsCMatrix) format for other population
#' @param sumstatsData summary data with required columns of SE_A, SE_B, A1.x, A1.y, and lFDR
#' @return returns a PRS specific LD matrix
#'
#' @import Matrix compiler
#'
#' @examples
#' sumstatsData = readRDS(file = system.file("extdata", "sumstatsData_toy.rds", package = "shaPRS") )
#'
#' # read SNP map files ( same toy data for the example)
#' pop1_map_rds = readRDS(file = system.file("extdata", "my_data.rds", package = "shaPRS") )
#' pop2_map_rds = readRDS(file = system.file("extdata", "my_data2.rds", package = "shaPRS") )
#'
#' # use chrom 21 as an example
#' chromNum=21
#'
#' # load the two chromosomes from each population ( same toy data for the example)
#' pop1LDmatrix = readRDS(file = system.file("extdata", "LDref.rds", package = "shaPRS") )
#' pop2LDmatrix = readRDS(file = system.file("extdata", "LDref2.rds", package = "shaPRS") )
#'
#'
#' # 2. grab the RSids from the map for the SNPS on this chrom,
#' # each LD mat has a potentiall different subset of SNPs
#' # this is guaranteed to be the same order as the pop1LDmatrix
#' pop1_chrom_SNPs = pop1_map_rds[ which(pop1_map_rds$chr == chromNum),]
#' # this is guaranteed to be the same order as the pop2LDmatrix
#' pop2_chrom_SNPs = pop2_map_rds[ which(pop2_map_rds$chr == chromNum),]
#' pop1_chrom_SNPs$pop1_id = 1:nrow(pop1_chrom_SNPs)
#' pop2_chrom_SNPs$pop2_id = 1:nrow(pop2_chrom_SNPs)
#'
#'
#' # intersect the 2 SNP lists so that we only use the ones common to both LD matrices by merging them
#' chrom_SNPs_df  <- merge(pop1_chrom_SNPs,pop2_chrom_SNPs, by = "rsid")
#'
#' # align the two LD matrices
#' chrom_SNPs_df = alignStrands(chrom_SNPs_df, A1.x ="a1.x", A2.x ="a0.x", A1.y ="a1.y", A2.y ="a0.y")
#'
#'
#' # align the summary for phe A and B
#' sumstatsData = alignStrands(sumstatsData)
#'
#' # subset sumstats data to the same chrom
#' sumstatsData = sumstatsData[which(sumstatsData$CHR == chromNum ),]
#'
#' # merge sumstats with common LD map data
#' sumstatsData  <- merge(chrom_SNPs_df,sumstatsData, by.x="rsid", by.y = "SNP")
#'
#' # remove duplicates
#' sumstatsData = sumstatsData[ !duplicated(sumstatsData$rsid) ,]
#' # use the effect alleles for the sumstats data with the effect allele of the LD mat
#' # as we are aligning the LD mats against each other, not against the summary stats
#' # we only use the lFDR /SE from the sumstats,
#' # which are directionless, so those dont need to be aligned
#' sumstatsData$A1.x =sumstatsData$a1.x
#' sumstatsData$A1.y =sumstatsData$a1.y
#'
#' # make sure the sumstats is ordered the same way as the LD matrix:
#' sumstatsData = sumstatsData[order(sumstatsData$pop1_id), ]
#' # it doesn't matter which matrix to use to order the sumstats as they are the same
#'
#' # subset the LD matrices to the SNPs we actualy have
#' pop1LDmatrix = pop1LDmatrix[sumstatsData$pop1_id,sumstatsData$pop1_id]
#' pop2LDmatrix = pop2LDmatrix[sumstatsData$pop2_id,sumstatsData$pop2_id]
#'
#' # generate the blended LD matrix
#' cormat = LDRefBlend(pop1LDmatrix,pop2LDmatrix, sumstatsData)
#'
#' # create a new map file that matches the SNPs common to both LD panels
#' map_rds_new = pop1_map_rds[which(pop1_map_rds$chr == chromNum),]
#' map_rds_new2 = map_rds_new[which(map_rds_new$rsid %in% sumstatsData$rsid),]
#'
#' # save the new LD matrix to a location of your choice
#' # saveRDS(cormat,file =paste0(<YOUR LOCATION>,"/LD_chr",chromNum,".rds"))
#'
#' # save its Map file too
#' #saveRDS(map_rds_new2,file = paste0(<YOUR LOCATION>,"/LD_chr",chromNum,"_map.rds"))
#'
#' @export
LDRefBlend = function(pop1LDmatrix,pop2LDmatrix, sumstatsData) {
  #library(compiler)
  enableJIT(3)

  #library(Matrix)

  ## use symbols from the tex
  wA = 1 - sumstatsData$lFDR    # wA is SNP_A's lFDR , between study 1 and study 2
  tauA1 = 1/sumstatsData$SE_A^2 # tauA1 is SNP_A's study 1 precision
  tauA2 = 1/sumstatsData$SE_B^2 # tauA2 is SNP_A's study 2 precision

  wB = 1 - sumstatsData$lFDR    # wB is SNP_B's lFDR, between study 1 and study 2
  tauB1 = 1/sumstatsData$SE_A^2 # tauB1 is SNP_B's study 1 precision
  tauB2 = 1/sumstatsData$SE_B^2 # tauB2 is SNP_B's study 2 precision

  ## build covariance matrix
  r1 = convertDSCToDense(pop1LDmatrix) # LD between SNPA and SNPB in study (population) 1 # must convert to dense otherwise will get memory error later
  prodA=(tauA1 + wA*tauA2)/(tauA1 + tauA2)/sqrt(tauA1) # vector of 5, IE pre-calculating all for all 'A' and all 'B; SNPs, all terms
  prodB=(tauB1 + wB*tauB2)/(tauB1 + tauB2)/sqrt(tauB1) # vector of 5
  # mat of 5x5, an outer product of A and B, IE make it the same dim as the LD matrix, IE the outer product expands out and performs all calculations in the loop in one go
  # IE outer product =~ loop
  term1=outer(prodA,prodB,"*") * r1   # multiplying elementwise by the LD
  remove(r1) # free up RAM

  # align the two LD matrices: invert correlations where one of the studies' effect allele is flipped wrt other studies' effect allele
  misalignedAlleleIndices = which( as.character(sumstatsData$A1.x) != as.character(sumstatsData$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
  if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned variants correlation were reversed" )
  # create a vector of the flipped alleles, 1 for aloigned, and -1 for misaligned
  flipped_mask =  rep(1, nrow(sumstatsData) )
  flipped_mask[misalignedAlleleIndices] = -1 # set
  # create a mask of the flipped alleles: this is just the outer product of the above vector
  flippedMat=outer(flipped_mask,flipped_mask,"*")  # create a matrix that can be used to flip correlations via elementwise multiplications


  r2 = convertDSCToDense(pop2LDmatrix) # LD between SNPA and SNPB in study (population) 2
  r2 = r2 * flippedMat # apply the flipping to the second LD mat
  remove(flippedMat) # free up RAM

  prodA=((1-wA) * tauA2)/(tauA1 + tauA2)/sqrt(tauA2)
  prodB=((1-wB) * tauB2)/(tauB1 + tauB2)/sqrt(tauB2)
  term2=outer(prodA,prodB,"*") * r2
  remove(r2) # free up RAM
  covmat=term1 + term2 # reconstruct full eq from 3.1

  # free up ram
  remove(term1)
  remove(term2)

  ## need variance to get cor
  # this replaces my
  # SEbar_A = sqrt( VarAbar )
  # VarAbar = wA^2 /taoA1 + (1-wA^2) / (taoA1 + taoA2)
  VbarA=wA^2/tauA1 + (1-wA^2)/(tauA1+tauA2) # vector of 5, of the variances
  VbarB=wB^2/tauB1 + (1-wB^2)/(tauB1+tauB2)
  outerProd = outer(sqrt(VbarA),sqrt(VbarB),"*") # this multiplies all combination of SEs together to create a 5x5 dim matrix
  # and thus replacing the loop that generates
  cormat=covmat / outerProd

  # free up ram
  remove(covmat)
  remove(outerProd)
  return(cormat)
}



#' Convert DSC to Dense matrix
#'
#' Helper function that converts an LD from sparse DSC to Dense format, in a given number of parts (to overcome RAM limitations)
#'
#' @param pop1LDmatrix LD reference matrix in RDS (dsCMatrix) format for target population
#' @param numparts (optional) how many parts should be used for converting the matrix
#' @return returns a dense LD matrix
#'
#' @export
convertDSCToDense = function(pop1LDmatrix, numparts = 3) {
  # convert dsc to dense matrix in a non-crashing way
  firstHalf = floor(ncol(pop1LDmatrix)/numparts)
  start = 0;
  end = 0
  r1 = NULL
  for (i in 1:numparts) {
    start = end
    if (i == numparts) { # last one has to go to the end
      end = ncol(pop1LDmatrix)
    } else {
      end = end + firstHalf
    }

    r1_p1 = as.matrix(pop1LDmatrix[,(start+1):end])
    r1 = cbind(r1,r1_p1)
    remove(r1_p1)
  }
  return(r1)
}
