##############################################################
# Blends summary statistics according to a continuous weighting scheme:
# it continuously blends the sub-phenotype and the combined phenotype summary statistics
# and generates an LDPred formatted 

# INPUT:
# (1) subphenoLoc: location of an LDPred formatted GWAS summary statistics file with the following signature:
#chr	pos	SNP	A1	A2	Freq1.Hapmap	b	se	p	N
#1	779322	rs4040617	G	A	X	-0.01553	0.008287	0.06099	17554
# (2) CombinedPhenoLoc: location of an LDPred formatted GWAS summary statistics file, with same signature
# (3) blendFactorLoc: a headerless 2 column file of the SNPid and lFDR, (produced by shaPRS_adjust.R)
# (4) outputLoc: location of the output files
# OUTPUT:
# (1) an LDPred formatted summary statistics file
##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")} 

subphenoLoc = args[1]
CombinedPhenoLoc = args[2]
blendFactorLoc = args[3]
outputLoc = args[4]



##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = T) 
CombinedPheno= read.table(CombinedPhenoLoc, header = T) 
blendingFactors= read.table(blendFactorLoc, header = F) 

# 2. merge the 3 tables together by RSid, so they are always aligned, x = subpheno  and    y = CombinedPheno
subpheno_CombinedPheno = merge(subpheno,CombinedPheno,by.x = "SNP",by.y = "SNP")
subpheno_CombinedPheno_blending = merge(subpheno_CombinedPheno,blendingFactors,by.x = "SNP",by.y = "V1")

# 3. Blend the subpheno and CombinedPheno together and create new summary statistics via following logic:
# blendingfactor = lFDR = prob that it is the null, that it was CombinedPheno NOT the subpheno
# basic formula: beta_phenoA x P(phenoA != phenoB) + beta_CombinedPheno x P(phenoA == phenoB)
# Var(aX+bY)=a^2 var(X) + b^2 var(Y)
# se is sqrt var, so:
# sqrt( Var(aX+bY) )= sqrt( a^2*var(X) + b^2*var(Y) )
blendedSE = sqrt(  (1-subpheno_CombinedPheno_blending$V2)^2 * subpheno_CombinedPheno_blending$se.x^2 +  subpheno_CombinedPheno_blending$se.y^2 * subpheno_CombinedPheno_blending$V2^2  )
blendedBeta=subpheno_CombinedPheno_blending$b.x * (1-subpheno_CombinedPheno_blending$V2) + subpheno_CombinedPheno_blending$b.y * subpheno_CombinedPheno_blending$V2
blendedp=2*pnorm( abs(blendedBeta)/blendedSE,lower=FALSE)

# 4. create new data frame to store the new summary stats
blendedSumstats = data.frame(subpheno_CombinedPheno_blending$chr.x, subpheno_CombinedPheno_blending$pos.x, subpheno_CombinedPheno_blending$SNP, subpheno_CombinedPheno_blending$A1.x, subpheno_CombinedPheno_blending$A2.x, subpheno_CombinedPheno_blending$Freq1.Hapmap.x,
                             blendedBeta,
                             blendedSE,
                             blendedp,
                             round(subpheno_CombinedPheno_blending$N.x * (1-subpheno_CombinedPheno_blending$V2) + subpheno_CombinedPheno_blending$N.y * subpheno_CombinedPheno_blending$V2)
                             )
colnames(blendedSumstats) = colnames(CombinedPheno)
blendedSumstats= blendedSumstats[match(subpheno$SNP, blendedSumstats$SNP),]
blendedSumstats <- na.omit(blendedSumstats) # remove any NAs of SNPs that couldn't be matched

# 5. write blended stats to disk
write.table(blendedSumstats, outputLoc, sep = "\t", row.names = F, col.names = T, quote = FALSE) 
print(paste("written blended sumstats to",outputLoc))


