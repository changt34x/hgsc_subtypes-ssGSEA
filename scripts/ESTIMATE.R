###############################################
# ESTIMATE purity testing for HGSC samples
#
# Chang, T
# ~~~~~~~~~~~~
# This script runs ESTIMATE based on existing
# Dataset and produces boxplots of gene
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
library(estimate)
library(ggplot2)
library(reshape2)

source('scripts/functions/LoadOVCA_Data.R')
source('scripts/functions/multiplot.R')

############################
# 2. Load Data
############################
# Load datasets specified in arguments
ExpData <- LoadOVCA_Data(datasets = args, genelist_subset = 'commongenes')

############################
# 3. Write Data for ESTIMATE
############################
outputDirectory <- 'output/dataset/'

for (dataset in 1:length(ExpData)) {
  write.table(ExpData[[dataset]], paste0(outputDirectory, names(ExpData)[dataset], '_processed.gct')
              , sep = '\t', col.names = NA)
}

############################
# 4. Compute ESTIMATE
############################
for (dataset in 1:length(ExpData)) {
  filterCommonGenes(input.f = paste0(outputDirectory, names(ExpData)[dataset], '_processed.gct'),
                    output.f = paste0('output/estimate/commongenes/', names(ExpData)[dataset], '_filtered.gct'),
                    id='GeneSymbol')
  if (dataset == 1 | dataset == 3 | dataset == 5) {
    estimateScore(paste0('output/estimate/commongenes/', names(ExpData)[dataset], '_filtered.gct'),
                  paste0('output/estimate/score/', names(ExpData)[dataset], '_estimate.gct'),
                  platform="affymetrix")
  } else {
    estimateScore(paste0('output/estimate/commongenes/', names(ExpData)[dataset], '_filtered.gct'),
                  paste0('output/estimate/score/', names(ExpData)[dataset], '_estimate.gct'),
                  platform="agilent")
  }
}

############################
# 4. Output histograms
############################
# Read ESTIMATE results
ESTIMATE_results <- list()

for (dataset in 1:length(ExpData)) {
  ESTIMATE_results[[dataset]] <- read.table(file = paste0('output/estimate/score/', names(ExpData)[dataset],
                                                          '_estimate.gct'),
             sep = '\t', skip = 2, header = T)
  rownames(ESTIMATE_results[[dataset]]) <- ESTIMATE_results[[dataset]][,2]
  ESTIMATE_results[[dataset]]$NAME <- NULL
  ESTIMATE_results[[dataset]]$Description <- NULL
  ESTIMATE_results[[dataset]] <- as.data.frame(t(ESTIMATE_results[[dataset]]))
}

# Plot as histogram
for (dataset in 1:length(ExpData)) {
  
  # Calculate Bins
  breaks <- pretty(range(ESTIMATE_results[[dataset]]$StromalScore), 
                   n = nclass.FD(ESTIMATE_results[[dataset]]$StromalScore), min.n = 1)
  bwidth <- breaks[2] - breaks[1]
  
  # Plot Stromal Score
  stromal <- ggplot(data = ESTIMATE_results[[dataset]], aes(x = StromalScore)) +
    geom_histogram(binwidth = bwidth, colour = "black", fill = "darkorange4") +
    labs(title = paste0("Stromal Score (", names(ExpData)[dataset], ")"),
         x = "Stromal Score", y = "Count")
  
  # Calculate Bins
  breaks <- pretty(range(ESTIMATE_results[[dataset]]$ImmuneScore), 
                   n = nclass.FD(ESTIMATE_results[[dataset]]$ImmuneScore), min.n = 1)
  bwidth <- breaks[2] - breaks[1]
  
  # Plot Immune Score
  immune <- ggplot(data = ESTIMATE_results[[dataset]], aes(x = ImmuneScore))+
    geom_histogram(binwidth = bwidth, colour = "black", fill = "deepskyblue4") +
    labs(title = paste0("Immune Score (", names(ExpData)[dataset], ")"),
        x = "Immune Score", y = "Count")
  
  # Calculate Bins
  breaks <- pretty(range(ESTIMATE_results[[dataset]]$ESTIMATEScore), 
                   n = nclass.FD(ESTIMATE_results[[dataset]]$ESTIMATEScore), min.n = 1)
  bwidth <- breaks[2] - breaks[1]
  
  # Plot ESTIMATE Score
  estimate <- ggplot(data = ESTIMATE_results[[dataset]], aes(x = ESTIMATEScore)) +
    geom_histogram(binwidth = bwidth, colour = "black", fill = "darkgoldenrod4") +
    labs(title = paste0("ESTIMATE Score (", names(ExpData)[dataset], ")"),
         x = "ESTIMATE Score", y = "Count")
  
  agg_hist <- multiplot(stromal, immune, estimate, cols = 1)
  
  # Save to PNG
  dev.copy(png, paste0('output/estimate/figures/', names(ExpData)[[dataset]], '.png'),
           width=8, height=10, units="in",res=500)
  agg_hist
  dev.off()
}

############################
# 5. Purity analysis
############################
# Analysis vectors
s_means <- c()
s_medians <- c()
s_sd <- c()
i_means <- c()
i_medians <- c()
i_sd <- c()
e_means <- c()
e_medians <- c()
e_sd <- c()
tp_mean <- c()
tp_median <- c()

# Calculate tumor purity for affymetrix
for (dataset in 1:length(ExpData)) {
  s_means <- c(s_means, mean(ESTIMATE_results[[dataset]]$StromalScore))
  s_medians <- c(s_medians, median(ESTIMATE_results[[dataset]]$StromalScore))
  s_sd <- c(s_sd, sd(ESTIMATE_results[[dataset]]$StromalScore))
  i_means <- c(i_means, mean(ESTIMATE_results[[dataset]]$ImmuneScore))
  i_medians <- c(i_medians, median(ESTIMATE_results[[dataset]]$ImmuneScore))
  i_sd <- c(i_sd, sd(ESTIMATE_results[[dataset]]$ImmuneScore))
  e_means <- c(e_means, mean(ESTIMATE_results[[dataset]]$ESTIMATEScore))
  e_medians <- c(e_medians, median(ESTIMATE_results[[dataset]]$ESTIMATEScore))
  e_sd <- c(e_sd, sd(ESTIMATE_results[[dataset]]$ESTIMATEScore))
  if (dataset == 1 | dataset == 3 | dataset == 5) {
    tp_mean <- c(tp_mean, mean(ESTIMATE_results[[dataset]]$TumorPurity))
    tp_median <- c(tp_median, median(ESTIMATE_results[[dataset]]$TumorPurity))
  } else {
    tp_mean <- c(tp_mean, NA)
    tp_median <- c(tp_median, NA)
  }
}

# Create dataframe of summary analysis
summary.df <- data.frame("Dataset" = names(ExpData),
                         "Stromal Mean" = s_means,
                         "Stromal Median" = s_medians,
                         "Stromal Std Dev" = s_sd,
                         "Immune Mean" = i_means,
                         "Immune Median" = i_medians,
                         "Immune Std Dev" = i_sd,
                         "ESTIMATE Mean" = e_means,
                         "ESTIMATE Median" = e_medians,
                         "ESTIMATE Std Dev" = e_sd,
                         "Tumor Purity Mean" = tp_mean,
                         "Tumor Purity Median" = tp_median)

# Graph summary
summary.gd <- melt(summary.df[,c(1,2,5,8)], id.vars = "Dataset")

ggplot(data = summary.gd, aes(x = Dataset, y = value)) +
  geom_bar(stat = "identity") +
  facet_grid(.~variable) +
  coord_flip() +
  labs(x = 'Score',y = 'Dataset')
