#!/bin/bash

exec &>hgsc_analysis.out 

############################################
# ssGSEA testing for ovarian cancer
# 
# Chang, T
# ~~~~~~~~~~~~~
# This script stores instructions to perform ssGSEA analysis
# of HGSC samples. To change the data sets that the analysis
# is performed on, change the arguments passed into Part 3.
# See the README for additional setup information.
# ~~~~~~~~~~~~~~~~~~~~~
############################################

################
# Part 0:
# Constants
################
DATA1="TCGA_eset"
DATA2="Mayo"
DATA3="GSE32062.GPL6480_eset"
DATA4="GSE9891_eset"
DATA5="GSE26712_eset"
DATA6=""
DATA7=""
DATA8=""

################
# Part 1:
# Install dependencies
################
Rscript --vanilla scripts/INSTALL.R

################
# Part 2:
# Extract LM22 signature
################
Rscript --vanilla scripts/extractLM22DEG.R

################
# Part 3:
# Run ssGSEA
################
R --no-save --args $DATA1 $DATA2 $DATA3 $DATA4 $DATA5 \
< scripts/ssGSEA.R

################
# Part 4:
# Run ESTIMATE
################
R --no-save --args $DATA1 $DATA2 $DATA3 $DATA4 $DATA5 \
< scripts/ESTIMATE.R