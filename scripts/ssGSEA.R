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
library(ggplot2)

source('scripts/functions/LoadOVCA_Data.R')
source('scripts/functions/ssGSEA_Aggregate_Table.R')

############################
# 2. Load Data
############################
# Load datasets specified in arguments
ExpData <- LoadOVCA_Data(datasets = args, genelist_subset = 'commongenes')

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
  
  # Create plots of ssGSEA scores by cell type and dataset
  for (i in 2:(dim(clustering)[2] + 1)) {
    for (j in (dim(clustering)[2] + 2):(dim(cluster_ssgsea)[2])) {
      png(filename = paste0(output_dir, '/', names(ExpData)[dataset], '_', data_label[j], ' k=',
                          i, '.png'))
      boxplot(cluster_ssgsea[,j]~cluster_ssgsea[,i], data = cluster_ssgsea, xlab = 'Subtypes',
              ylab = 'Abundance', main = paste0(names(ExpData)[dataset], ' ', data_label[j],
                                                ' k = ', i))
      dev.off()
    }
  }
}

############################
# 4. Aggregate analysis (Cell Type x Dataset)
############################
# Get aggregate dataframe
agg_data <- ssGSEA_Aggregate_Table(args, dim(clustering)[2])

# Get types of cells analyzed
cell_types <- unique(agg_data$Cell.Type)

# Create output directory
dir.create('output/figures/lm22_overall', showWarnings = F)

# Plot data
for (i in 2:(dim(clustering)[2] + 1)) {
  # Create output directory by subtyping number
  output_dir <- paste0('output/figures/lm22_overall/k=', i)
  dir.create(output_dir, showWarnings = F)
  
  for (j in 1:length(cell_types)) {
    # Get cell type
    current_type <- cell_types[j]
    
    # Select data from table
    sel <- agg_data$K == i & agg_data$Cell.Type == current_type
    
    # Plot and save as
    ggplot(data = agg_data[sel,], mapping = aes(x = as.character(agg_data$Subtype[sel]), 
                                                y = agg_data$ssGSEA[sel], 
                                                fill = agg_data$Dataset[sel])) +
      geom_boxplot() +
      labs(title = paste0("lm22: ", current_type, " (K = ", i, ")"), x = "Subtype", 
           y = "Abundance") +
      guides(fill = guide_legend(title="Dataset"))
    
    ggsave(filename = paste0(output_dir, '/', current_type, '_k=', i, '.png'),
           plot = last_plot())
  }
}

############################
# 5. Aggregate analysis (Dataset x Cell Type)
############################
# Plot data
for (dataset in 1:length(ExpData)) {
  dataset_name <- names(ExpData)[dataset]
  
  output_dir <- paste0('output/figures/', dataset_name, '/summary')
  dir.create(output_dir, showWarnings = F)
  
  for (i in 2:(dim(clustering)[2] + 1)) {
    for (j in 1:i) {
      # Select data from table
      sel <- agg_data$K == i & agg_data$Dataset == dataset_name & agg_data$Subtype == j
      
      # Plot and save as png
      ggplot(data = agg_data[sel,], mapping = aes(x = agg_data$Cell.Type[sel],
                                                  y = agg_data$ssGSEA[sel])) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = paste0(dataset_name, " (K = ", i, ", Subtype = ", j, ")"),
             x = "Cell Type", y = "Abundance")
      
      ggsave(filename = paste0(output_dir, '/', dataset_name, '_k=', i, '_s=', j, '.png'),
             plot = last_plot())
    }
  }
}

for (i in 2:(dim(clustering)[2] + 1)) {
  # Create output directory by subtyping number
  output_dir <- paste0('output/figures/lm22_overall/k=', i)
  dir.create(output_dir, showWarnings = F)
  
  for (j in 1:length(cell_types)) {
    # Get cell type
    current_type <- cell_types[j]
    
    # Select data from table
    sel <- agg_data$K == i & agg_data$Cell.Type == current_type
    
    # Plot and save as png
    ggplot(data = agg_data[sel,], mapping = aes(x = as.character(agg_data$Subtype[sel]), 
                                                y = agg_data$ssGSEA[sel], 
                                                fill = agg_data$Dataset[sel])) +
      geom_boxplot() +
      labs(title = paste0("lm22: ", current_type, " (K = ", i, ")"), x = "Subtype", 
           y = "Abundance") +
      guides(fill = guide_legend(title="Dataset"))
    
    ggsave(filename = paste0(output_dir, '/', current_type, '_k=', i, '.png'),
           plot = last_plot())
  }
}