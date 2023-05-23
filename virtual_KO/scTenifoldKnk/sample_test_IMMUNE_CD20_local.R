# libraries ---------------------------------------------------------------
# install.packages("scTenifoldKnk")
library(scTenifoldKnk)
library(tidyverse)

# read in the data --------------------------------------------------------
# this is a sample dataset from the 
# scRNAseq_loc <- system.file("single-cell/example.csv",package="scTenifoldKnk")
# scRNAseq <- read.csv(scRNAseq_loc, row.names = 1)
scRNAseq <- read_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_raw_data_all.tsv") %>%
  column_to_rownames(var = "gene")

# explore the dataset
dim(scRNAseq)

# the columns are the barcodes and the rows are the genes
# the table is a raw table of counts
scRNAseq[1:10,1:10]

# verify the presence of the genes of interest
rownames(scRNAseq) %>%
  str_subset(pattern = "MS4A1")

# verify the expression of the two markers genes in the dataset
rowSums(scRNAseq)[c("MS4A1")]
# number of positive cells
rowSums(scRNAseq>0)[c("MS4A1")]

# run the tool ------------------------------------------------------------
# one gene
test_CD20 <- scTenifoldKnk(countMatrix = scRNAseq, gKO = c("CD20"), qc_minLSize = 0)
saveRDS(test_CD20,"out/IMMUNE_CCA_martina_trimm_KO_CD20.rds")