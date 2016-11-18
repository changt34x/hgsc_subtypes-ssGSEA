###############################################
# ssGSEA testing for ovarian cancer
#
# Chang, T
# ~~~~~~~~~~~~
# This script documents the install process for
# required tools.
###############################################

# Install standard R packages
mirror <- "http://cran.us.r-project.org"
packages <- c("rJava","xlsx", "estimate")
install.packages(packages, repos = mirror, dependencies = T)

# Install Estimate
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

###############################################
# Install BioConductor Tools
###############################################
source("https://bioconductor.org/biocLite.R")
bioc_packages <- c("GSEABase", "GSVA", "curatedOvarianData")
biocLite(bioc_packages)