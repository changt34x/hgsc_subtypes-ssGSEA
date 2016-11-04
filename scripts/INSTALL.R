###############################################
# ssGSEA testing for ovarian cancer
#
# Chang, T
# ~~~~~~~~~~~~
# This script documents the install process for
# required tools.
###############################################

mirror <- "http://cran.us.r-project.org"
packages <- c("rJava","xlsx")
install.packages(packages, repos = mirror)

###############################################
# Install BioConductor Tools
###############################################
source("https://bioconductor.org/biocLite.R")
bioc_packages <- c("GSEABase", "GSVA", "curatedOvarianData")
biocLite(bioc_packages)
