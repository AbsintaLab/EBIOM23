# libraries ---------------------------------------------------------------
library(scTenifoldNet)
library(tidyverse)

# read in the data --------------------------------------------------------
# raw matrix all the cells
IMM_ALL <- read_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_ALL.tsv") %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(IMM_ALL)
IMM_ALL[1:10,1:10]

# raw matrix removed the Bcells
IMM_noCD20cell <- read_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeCD20cells.tsv") %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(IMM_noCD20cell)
IMM_noCD20cell[1:10,1:10]

# run the method ----------------------------------------------------------
output <- scTenifoldNet(X = IMM_ALL, Y = IMM_noCD20cell)

# save the output ---------------------------------------------------------
saveRDS(output,"out/object/IMMUNE_CCA_martina_scTenifoldNet_CD20cells.rds")
