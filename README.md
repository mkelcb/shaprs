# shaprs
ShaPRS: Leveraging shared genetic effects across traits and ancestries improves accuracy of polygenic scores

Installation:

``` R
install_github("mkelcb/shaprs")
library("shaPRS")
```

## To find the shaPRS weighted meta-analysis of a proximal and adjunct data, simply run:
``` R
proximalLoc <- system.file("extdata", "phenoA_sumstats", package = "shaPRS")
adjunctLoc <- system.file("extdata", "phenoB_sumstats", package = "shaPRS")
shaPRS(proximalLoc, adjunctLoc, "<YOUR_OUTPUT_FOLDER>")
``` 
- This will output your final summary statistics file with the postfix "_shaprs" that you may use in your favourite PRS generation tool. 
- The above will also output a few other files that may be of interest: "_meta"  (fixed fixed effect meta-analysis) and "_SNP_lFDR" (lFDR estimates and Q-values for each SNP). 


## Blend LD ref matrices (optional):

``` R
sumstatsData = readRDS(file = system.file("extdata", "sumstatsData_toy.rds", package = "shaPRS") )
```


read SNP map files ( same toy data for the example)
``` R
pop1_map_rds = readRDS(file = system.file("extdata", "my_data.rds", package = "shaPRS") )
pop2_map_rds = readRDS(file = system.file("extdata", "my_data2.rds", package = "shaPRS") )
```

use chrom 21 as an example
``` R
chromNum=21
```

load the two chromosomes from each population ( same toy data for the example)

``` R
pop1LDmatrix = readRDS(file = system.file("extdata", "LDref.rds", package = "shaPRS") )
pop2LDmatrix = readRDS(file = system.file("extdata", "LDref2.rds", package = "shaPRS") )
```


2. grab the RSids from the map for the SNPS on this chrom,
each LD mat has a potentiall different subset of SNPs
this is guaranteed to be the same order as the pop1LDmatrix

```
pop1_chrom_SNPs = pop1_map_rds[ which(pop1_map_rds$chr == chromNum),]
```

this is guaranteed to be the same order as the pop2LDmatrix

```
pop2_chrom_SNPs = pop2_map_rds[ which(pop2_map_rds$chr == chromNum),]
pop1_chrom_SNPs$pop1_id = 1:nrow(pop1_chrom_SNPs)
pop2_chrom_SNPs$pop2_id = 1:nrow(pop2_chrom_SNPs)
```


intersect the 2 SNP lists so that we only use the ones common to both LD matrices by merging them
```
chrom_SNPs_df  <- merge(pop1_chrom_SNPs,pop2_chrom_SNPs, by = "rsid")
```

align the two LD matrices
```
chrom_SNPs_df = alignStrands(chrom_SNPs_df, A1.x ="a1.x", A2.x ="a0.x", A1.y ="a1.y", A2.y ="a0.y")
```


align the summary for phe A and B
```
sumstatsData = alignStrands(sumstatsData)
```

subset sumstats data to the same chrom
```
sumstatsData = sumstatsData[which(sumstatsData$CHR == chromNum ),]
```

merge sumstats with common LD map data
```
sumstatsData  <- merge(chrom_SNPs_df,sumstatsData, by.x="rsid", by.y = "SNP")
```

remove duplicates
```
sumstatsData = sumstatsData[ !duplicated(sumstatsData$rsid) ,]
```

use the effect alleles for the sumstats data with the effect allele of the LD mat
as we are aligning the LD mats against each other, not against the summary stats
we only use the lFDR /SE from the sumstats,
which are directionless, so those dont need to be aligned

```
sumstatsData$A1.x =sumstatsData$a1.x
sumstatsData$A1.y =sumstatsData$a1.y
```

make sure the sumstats is ordered the same way as the LD matrix:
```
sumstatsData = sumstatsData[order(sumstatsData$pop1_id), ]
```

(it doesn't matter which matrix to use to order the sumstats as they are the same)

subset the LD matrices to the SNPs we actualy have
```
pop1LDmatrix = pop1LDmatrix[sumstatsData$pop1_id,sumstatsData$pop1_id]
pop2LDmatrix = pop2LDmatrix[sumstatsData$pop2_id,sumstatsData$pop2_id]
```

generate the blended LD matrix
```
cormat = LDRefBlend(pop1LDmatrix,pop2LDmatrix, sumstatsData)
```

create a new map file that matches the SNPs common to both LD panels
>map_rds_new = pop1_map_rds[which(pop1_map_rds$chr == chromNum),]

>map_rds_new2 = map_rds_new[which(map_rds_new$rsid %in% sumstatsData$rsid),] 

save the new LD matrix to a location of your choice

>saveRDS(cormat,file =paste0(\<YOUR LOCATION\>,"/LD_chr",chromNum,".rds"))

save its Map file too

>saveRDS(map_rds_new2,file = paste0(\<YOUR LOCATION\>,"/LD_chr",chromNum,"_map.rds"))

- The cormat is a 29x29 dense matrix of SNP-SNP correlations, which are saved to a location of your choice, together with its map file.
