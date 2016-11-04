# High-Grade Serous Ovarian Cancer Subtypes - Cellular Deconvolution Through ssGSEA

## Summary
This repository contains a procedure that can be used to apply ssGSEA analysis to 
existing gene expression and subtyping data for hgsc datasets. The experimental
data set used in this example was produced by The Cancer Genome Atlas (2011) and
is available from BioConductor's curatedOvarianData package. The same procedure
should be applicable and easily adaptible to other hgsc datasets and those from
other diseases.

## Reproduction

### Step 1: Download signature gene file
You will need Supplementary Table S1 from
[Newman et al. 2015](https://doi.org/10.1038/nmeth.3337).

Please download the file. Copy A3:W550 and paste into A1 of a new spreadsheet.
Save this file as '/data/raw/LM22_DEG.csv'.

### Step 2: Run analysis pipeline
Run ANALYSIS.sh