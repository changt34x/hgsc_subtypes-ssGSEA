###############################################
# ssGSEA testing for ovarian cancer
#
# Chang, T
# ~~~~~~~~~~~~
# This script runs ssGSEA based on existing
# subtype data and produces boxplots of gene
# abundances
###############################################

############################
# 0. Parse Arguments
############################
args <- commandArgs(trailingOnly = T)
#args <- c("TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset")

############################
# 1. Load Libraries
############################
library(GSVA)
library(curatedOvarianData)

source('scripts/functions/LoadOVCA_Data.R')

############################
# 2. Load Data
############################
# Load datasets specified in arguments
ExpData <- LoadOVCA_Data(datasets = args, genelist_subset = 'None')

# Load lm22 gene set
lm22_genes <- read.table('data/lm22_genes.tsv', sep = '\t', header = T, row.names = 1)
lm22_geneset <- lapply(lm22_genes, function(x) rownames(lm22_genes)[x == 1])

############################
# 3. Run ssGSEA
############################
for (dataset in 1:length(ExpData)) {
  
  # Load expression data
  gene_exp <- ExpData[[dataset]]
  gene_input <- rownames(gene_exp)

  # Load clustering data
  clustering <- read.csv(paste0('data/ClusterMembership/', names(ExpData)[dataset], '_cluster.csv')
                         , stringsAsFactors = F, row.names = 1)
  
  # Run ssGSEA function
  ssgsea_result <- GSVA::gsva(expr = gene_exp, gset.idx.list = lm22_geneset,
                              method = 'ssgsea', verbose = FALSE)

  # Combine ssGSEA data and clustering
  cluster_ssgsea <- merge(clustering, t(ssgsea_result), by=c('row.names'), all = F)

  # Write results
  write.table(cluster_ssgsea, paste0('output/tables/', names(ExpData)[dataset], '_ssGSEA.tsv')
              , sep = '\t', col.names = NA)

  data_label <- colnames(cluster_ssgsea)

  output_dir <- paste0('output/figures/', names(ExpData)[dataset])
  dir.create(output_dir, showWarnings = F)
  
  for (i in 2:4) {
    for (j in 5:(dim(cluster_ssgsea)[2])) {
      png(filename=paste0(output_dir, '/', names(ExpData)[dataset], '_', data_label[j], ' k=',
                          i, '.png'))
      boxplot(cluster_ssgsea[,j]~cluster_ssgsea[,i], data = cluster_ssgsea, xlab = 'Subtypes',
              ylab = 'Abundance', main = paste0(names(ExpData)[dataset], ' ', data_label[j],
                                                ' k = ', i))
      dev.off()
    }
  }
}
