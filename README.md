# shaprs
Leveraging shared genetic effects to improve  genetic risk prediction for related diseases

The process relies on two simple, stand alone R scripts that requires no installation

## Step 1: adjust your summary statistics, run:
>shaPRS_adjust.R shapersToydata.txt 0 myOutput 0.5 0.99

this will produce 3 files:

**(1) myOutput_SNP_lFDR**: a table of lFDRs estimate for each SNP's heterogeneity

**(2) myOutput_SNPs_0.5**: SNPs that are found to be heterogeneous at an lFDR threshold of 0.5

**(3) myOutput_SNPs_0.99**: SNPs that are found to be heterogeneous at an lFDR threshold of 0.99

- if your goal is to simply create a Composite PRS, take either files (2-3) as a list of SNPs for which the subphenotype estimates  should be used

- if you want more thresholds, you can add as many as you would like


## Step 2: Blend summary statistics according to a continuous weighting scheme, run:
>sumstatsBlender_shaPRS.R phenoA_sumstats Combined_sumstats myOutput_SNP_lFDR myOutput

this will produce 1 file:

**(1) myOutput**: LDPred formatted summary statistics file

That's it. You may now then use the new file in LDPred or the favourite PRS generation tool. 

If you are interested, take a look at the .R scripts, which have more comments in them.
