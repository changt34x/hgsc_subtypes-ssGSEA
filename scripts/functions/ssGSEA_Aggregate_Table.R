ssGSEA_Aggregate_Table <- function(datasets, clusters = 3) {
  
  # Create summary vectors
  SampleID <- c()
  CellType <- c()
  Subtype <- c()
  Dataset <- c()
  K <- c()
  Score <- c()
  
  for (i in 1:length(datasets)) {
    
    # Read dataset ssGSEA analysis
    analysis_data <- read.table(paste0('output/tables/', datasets[i], '_ssGSEA.tsv'),
                               sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
    
    # Get number of samples
    row_count <- dim(analysis_data)[1]
    
    # Run for each k clustering and cell type
    for (j in 2:(clusters + 1)) {
      for (k in (clusters + 2):dim(analysis_data)[2]) {
        
        # Add data to summary vectors
        SampleID <- c(SampleID, analysis_data[,1])
        CellType <- c(CellType, rep(colnames(analysis_data)[k], row_count))
        Subtype <- c(Subtype, analysis_data[, j])
        Dataset <- c(Dataset, rep(datasets[i], row_count))
        K <- c(K, rep(max(analysis_data[, j]), row_count))
        Score <- c(Score, analysis_data[, k])
      }
    }
  }
  
  # Create summary dataframe
  summary_data <- data.frame("SampleID" = SampleID, "Cell Type" = CellType, 
                             "Subtype" = Subtype, "Dataset" = Dataset,
                             "K" = K, "ssGSEA" = Score)
  
  # Return summary dataframe
  return(summary_data)
}