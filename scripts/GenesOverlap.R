library(stringr)
library(stringi)
library(gdata)

DATA1 <- "TCGA_HumanMethylation27"
DATA2 <- "LM22.txt"
GENES <- "CommonGenes_genelist.csv"

' Exp1 <- read.table(paste("data/raw/", DATA1, sep=""), header = T, stringsAsFactors = T) '
Exp2 <- read.table(paste("data/raw/", DATA2, sep=""), header = F, stringsAsFactors = F)
GeneTable <- read.table(paste("data/raw/", GENES, sep=""), header = T, stringsAsFactors = F)

' Exp1Gene <- Exp1[[1]] '
Exp2Gene <- Exp2[[1]]
GoodGenes <- GeneTable[[1]]

# Conform strings from Exp data set
trim.conform <- function(x) {
  temp <- substring(x, 5)
  temp <- str_extract(temp, "\\w{3}-\\w+")
  return(toupper(temp))
}

trim.dash <- function(x) {
  part1 <- substring(x, 1, 3)
  part2 <- substring(x, 5)
  return(paste(part1, part2, sep=""))
}

Exp2Gene <- trim.conform(Exp2Gene)
Exp2Gene <- trim.dash(Exp2Gene)

# Conform strings from good genes
GoodGenesFilter <- c()

for (i in 1:length(GoodGenes)) {
  TempGene <- GoodGenes[i]
  toRun = TRUE
  while (toRun) {
    if (grepl("\\w+///", TempGene) == TRUE) {
      TempGeneSec <- str_extract(TempGene, "\\w+///")
      toExclude <- stri_length(TempGeneSec)
      TempGeneSec <- substring(TempGeneSec, 1, toExclude  - 3)
      TempGene <- substring(TempGene, toExclude + 1)
      GoodGenesFilter <- c(GoodGenesFilter, TempGeneSec)
    } else {
      toRun = FALSE
    }
  }
  GoodGenesFilter <- c(GoodGenesFilter, TempGene)
}

# Match genes with good genes
Match <- c()

for (i in 1:length(Exp2Gene)) {
  if (Exp2Gene[i] %in% GoodGenesFilter) {
    Match <- c(Match, Exp2Gene[i])
  }
}

# 446 Genes overlapping with Good Genes