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
IMM_noTcell <- read_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeTcells.tsv") %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(IMM_noTcell)
IMM_noTcell[1:10,1:10]

# run the method ----------------------------------------------------------
output <- scTenifoldNet(X = IMM_ALL, Y = IMM_noTcell)

# save the output ---------------------------------------------------------
saveRDS(output,"out/object/IMMUNE_CCA_martina_scTenifoldNet_Tcells.rds")
