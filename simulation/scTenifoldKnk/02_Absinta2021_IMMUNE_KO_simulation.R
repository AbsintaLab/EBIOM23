# libraries ---------------------------------------------------------------
library(scTenifoldKnk)
library(tidyverse)

# read in the data --------------------------------------------------------
scRNAseq <- read_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_raw_data_all.tsv") %>% 
  column_to_rownames(var = "gene")

# explore the dataset
dim(scRNAseq)

# the columns are the barcodes and the rows are the genes
# the table is a raw table of counts
scRNAseq[1:10,1:10]

# verify the presence of the genes of interest
rownames(scRNAseq) %>% str_subset(pattern = "BTK")

# verify the expression of the two markers genes in the dataset
rowSums(scRNAseq)[c("BTK")]
# number of positive cells
rowSums(scRNAseq>0)[c("BTK")]

# run the tool ------------------------------------------------------------
# one gene
test_BTK <- scTenifoldKnk(countMatrix = scRNAseq, gKO = c("BTK"), qc_minLSize = 0)
saveRDS(test_BTK,"out/IMMUNE_CCA_martina_trimm_KO_BTK.rds")
