###############################################
# ssGSEA testing for ovarian cancer
#
# Chang, T
# ~~~~~~~~~~~~
# This script takes Supplmentary Table S1 from
# Newman et al. 2015 and prepares it into a csv
# format while removing common genes to be
# tested.
###############################################

# Read raw LM22 data
lm22_genes <- read.csv('data/raw/LM22_DEG.csv', header = T, row.names = 1)

# Read common genes data
common_genes <- read.csv('data/Genes/CommonGenes_genelist.csv', header = T, row.names = 1)

# Gene selection
lm22Row <- rownames(lm22_genes)
commonRow <- rownames(common_genes)
similar <- lm22Row[lm22Row %in% commonRow]
lm22_filtered <- lm22_genes[similar, ]

# Output table
write.table(lm22_filtered, 'data/lm22_genes.tsv', sep = '\t', col.names = NA)
