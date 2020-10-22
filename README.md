# shaprs
Leveraging shared genetic effects to improve  genetic risk prediction for related diseases

Installation:
>install_github("mkelcb/shaprs")
>library("shaPRS")
## Step 1: adjust your summary statistics, run:
>inputDataLoc <- system.file("extdata", "shapersToydata.txt", package = "shaPRS")
>inputData= read.table(inputDataLoc, header = T)
>results = shaPRS_adjust(inputData, thresholds=c(0.5,0.99))

- the results object will have two lists,'lFDRTable', which provides the lFDR estimates for each SNP and 'hardThresholds' which itself is a list of the SNPs that found to be heterogeneous at the thresholds 0.5 and 0.99. 

- if your goal is to simply create a Composite PRS, take either files (2-3) as a list of SNPs for which the subphenotype estimates  should be used

- if you want more thresholds, you can add as many as you would like


## Step 2: Blend summary statistics according to a continuous weighting scheme, run:
>subphenoLoc <- system.file("extdata", "phenoA_sumstats", package = "shaPRS")
>CombinedPhenoLoc <- system.file("extdata", "Combined_sumstats", package = "shaPRS")
>blendFactorLoc <- system.file("extdata", "myOutput_SNP_lFDR", package = "shaPRS")
>subpheno= read.table(subphenoLoc, header = T)
>CombinedPheno= read.table(CombinedPhenoLoc, header = T)
>blendingFactors= read.table(blendFactorLoc, header = F)
>blendedSumstats = shaPRS_blend(subpheno, CombinedPheno, blendingFactors)

- 'blendedSumstats' is an LDPred formatted summary statistics file

That's it. You may now then use the new file in LDPred or the favourite PRS generation tool. 
