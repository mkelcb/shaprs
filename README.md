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


## Blend LD ref matrices (cross-ancestry analysis):

``` R
Pop1LDRefLoc <- paste0(system.file("extdata", "", package = "shaPRS"), "/")
Pop2LDRefLoc <- paste0(system.file("extdata", "", package = "shaPRS"), "/")
blendFactorLoc <- system.file("extdata", "pop_SNP_lFDR", package = "shaPRS")
adjustinputLoc <- system.file("extdata", "pop_adjustinput", package = "shaPRS")
outputLoc <- "<YOUR LOCATION>"
shaPRS_LDGen(Pop1LDRefLoc, Pop2LDRefLoc, blendFactorLoc, adjustinputLoc, outputLoc)
```

- This runs the shaPRS cross-ancestry analysis on the included toy dataset and generates the LD-reference panel for the 22 autosomes, along with a map.rds file.
- To run it on real data, first run the main shaPRS() to generate the "_SNP_lFDR" and "_adjustinput" files (see above) for "blendFactorLoc" and "adjustinputLoc". Finally, specify two LDpred2 formatted LD-reference panels, appropriate to your proximal/adjunct datasets ("Pop1LDRefLoc" and "Pop2LDRefLoc").
- The output data in the results folder (<YOUR LOCATION>), can then be used by LDpred2, processed by the same scripts that you have been using on the standard LD-ref matrices.
