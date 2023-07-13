# libraries ---------------------------------------------------------------
library(tidyverse)

# read in the data --------------------------------------------------------
meta <- read_tsv("../scTenifoldKnk/data/IMMUNE_martina_meta.tsv")

# define the barcodes of the plasmacells
id_Bcells <- meta %>%
  filter(seurat_clusters == 10,
         pathology == "c_chronic_active") %>%
  pull(barcode)

length(id_Bcells)

meta %>%
  filter(barcode %in% id_Bcells)

# define the barcodes of the CD4 cells
id_Tcells <- meta %>%
  filter(seurat_clusters == 7,
         pathology == "c_chronic_active") %>%
  pull(barcode)

length(id_Tcells)

meta %>%
  filter(barcode %in% id_Tcells)

# define the barcodes of cluster 1 (mims)
id_mim1 <- meta %>%
  filter(seurat_clusters == 1,
         pathology == "c_chronic_active") %>%
  pull(barcode)

length(id_mim1)

meta %>%
  filter(barcode %in% id_mim1)

# define the barcodes of cluster 8 (mims)
id_mim8 <- meta %>%
  filter(seurat_clusters == 8,
         pathology == "c_chronic_active") %>%
  pull(barcode)

length(id_mim8)

meta %>%
  filter(barcode %in% id_mim8)

# load the expression data from CCA and remove all the gene expressed in less than 5% of the cells
df_IMM_CCA <- read_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_raw_data_all.tsv")
df_IMM_CCA_norm <- read_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_norm_data_all.tsv")
df_ALL_CCA <- read_tsv("../out_large/scTenifoldKnk/table/ALL_CCA_martina_raw_data_all.tsv")

dim(df_IMM_CCA)
dim(df_IMM_CCA_norm)
dim(df_ALL_CCA)

# how many genes are trimmed
sum((rowSums(df_IMM_CCA[,-1]>0)/dim(df_IMM_CCA)[2])>0.05)
sum((rowSums(df_IMM_CCA_norm[,-1]>0)/dim(df_IMM_CCA_norm)[2])>0.05)
sum((rowSums(df_ALL_CCA[,-1]>0)/dim(df_ALL_CCA)[2])>0.05)

id_genes_IMM_CCA <- (rowSums(df_IMM_CCA[,-1]>0)/dim(df_IMM_CCA)[2])>0.05
id_genes_IMM_CCA_norm <- (rowSums(df_IMM_CCA_norm[,-1]>0)/dim(df_IMM_CCA_norm)[2])>0.05
id_genes_ALL_CCA <- (rowSums(df_ALL_CCA[,-1]>0)/dim(df_ALL_CCA)[2])>0.05

# check the expression of CD20 and simulate the removal of the only those non zero expression cells
df_IMM_CCA$gene %>%
  str_subset(pattern = "MS4A1")

# verify the expression of the two markers genes in the dataset
rowSums(df_IMM_CCA[,-1])[df_IMM_CCA$gene == "MS4A1"]
# number of positive cells
rowSums(df_IMM_CCA[,-1]>0)[df_IMM_CCA$gene == "MS4A1"]
# pull the barcodes of the 9 positive cells
id_CD20cells <- df_IMM_CCA[df_IMM_CCA$gene == "MS4A1",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(exp>0) %>%
  pull(barcodes)

id_CD20cells_norm <- df_IMM_CCA_norm[df_IMM_CCA_norm$gene == "MS4A1",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(exp>0) %>%
  pull(barcodes)

# test if other genes are coexpressed
test <- df_IMM_CCA[,c("gene",id_CD20cells)] %>% 
  filter(gene %in% c("CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C"))

rowSums(test[,-1]>0)

meta %>%
  filter(barcode %in% id_CD20cells)

meta %>%
  filter(barcode %in% id_CD20cells_norm)

df_IMM_CCA[df_IMM_CCA$gene == "MS4A1",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(barcodes %in% id_CD20cells)

df_IMM_CCA_norm[df_IMM_CCA_norm$gene == "MS4A1",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(barcodes %in% id_CD20cells)

# check the expression of BTK
df_IMM_CCA$gene %>%
  str_subset(pattern = "BTK")

# verify the expression of the two markers genes in the dataset
rowSums(df_IMM_CCA[,-1])[df_IMM_CCA$gene == "BTK"]
# number of positive cells
rowSums(df_IMM_CCA[,-1]>0)[df_IMM_CCA$gene == "BTK"]
# pull the barcodes of the 9 positive cells
id_BTKcells <- df_IMM_CCA[df_IMM_CCA$gene == "BTK",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(exp>0) %>%
  pull(barcodes)

meta %>%
  filter(barcode %in% id_BTKcells)

summ_all <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(tot_cluster = n())

meta %>%
  filter(barcode %in% id_BTKcells) %>%
  group_by(seurat_clusters) %>%
  summarise(BTK_pos = n()) %>%
  left_join(summ_all,by = "seurat_clusters") %>%
  mutate(prop = BTK_pos/tot_cluster)

df_IMM_CCA[df_IMM_CCA$gene == "BTK",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(barcodes %in% id_BTKcells)

# check the removal of CD69 from T cells
# the gene is in the dataset
df_IMM_CCA$gene %>%
  str_subset(pattern = "CD69")

# verify the expression of the two markers genes in the dataset
rowSums(df_IMM_CCA[,-1])[df_IMM_CCA$gene == "CD69"]
# number of positive cells
rowSums(df_IMM_CCA[,-1]>0)[df_IMM_CCA$gene == "CD69"]

# pull the barcodes of the positive cells
id_CD69cells <- df_IMM_CCA[df_IMM_CCA$gene == "CD69",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(exp>0) %>%
  pull(barcodes)

id_CD69cells_norm <- df_IMM_CCA_norm[df_IMM_CCA_norm$gene == "CD69",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(exp>0) %>%
  pull(barcodes)

# enumerating the CD69 positive
meta %>%
  filter(barcode %in% id_CD69cells)
# filter the CD69 positive and T celkls
meta %>%
  filter(barcode %in% id_CD69cells) %>%
  filter(barcode %in% id_Tcells)

meta %>%
  filter(barcode %in% id_CD69cells_norm) %>%
  filter(barcode %in% id_Tcells)

df_IMM_CCA[df_IMM_CCA$gene == "CD69",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(barcodes %in% id_CD69cells) %>%
  filter(barcodes %in% id_Tcells)

df_IMM_CCA_norm[df_IMM_CCA_norm$gene == "CD69",] %>%
  pivot_longer(names_to = "barcodes",values_to = "exp",-gene) %>%
  filter(barcodes %in% id_CD69cells) %>%
  filter(barcodes %in% id_Tcells)

# -------------------------------------------------------------------------
# save the trimmed version of the matrices after trimming the genes and the cells
# IMM
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_ALL.tsv")

# save the normalized version of the tables
df_IMM_CCA_norm %>%
  filter(id_genes_IMM_CCA_norm) %>%
  dim()

df_IMM_CCA_norm %>%
  filter(id_genes_IMM_CCA_norm) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_ALL_norm.tsv")

# remove the B cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_Bcells)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_Bcells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeBcells.tsv")

# remove the mim1 cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_mim1)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_mim1)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeMIM1cells.tsv")

# remove the mim8 cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_mim8)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_mim8)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeMIM8cells.tsv")

# remove the CD20 cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_CD20cells)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_CD20cells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeCD20cells.tsv")

# remove the BTK cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  dplyr::select(-all_of(id_BTKcells)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  dplyr::select(-all_of(id_BTKcells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeBTKcells.tsv")

# remove the B cells non CD20 cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(setdiff(id_Bcells,id_CD20cells))) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_CD20cells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeCD20cells.tsv")

# remove the T cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_Tcells)) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(id_Tcells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeTcells.tsv")

# ALL
df_ALL_CCA %>%
  filter(id_genes_ALL_CCA) %>%
  dim()

df_ALL_CCA %>%
  filter(id_genes_ALL_CCA) %>%
  write_tsv("../out_large/scTenifoldNet/table/ALL_CCA_martina_Net_trimm_ALL.tsv")

# remove the B cells
df_ALL_CCA %>%
  filter(id_genes_ALL_CCA) %>%
  select(-all_of(id_Bcells)) %>%
  dim()

df_ALL_CCA %>%
  filter(id_genes_ALL_CCA) %>%
  select(-all_of(id_Bcells)) %>%
  write_tsv("../out_large/scTenifoldNet/table/ALL_CCA_martina_Net_trimm_removeBcells.tsv")

# remove the T cells CD69 positive cells
df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(intersect(id_Tcells,id_CD69cells))) %>%
  dim()

df_IMM_CCA %>%
  filter(id_genes_IMM_CCA) %>%
  select(-all_of(intersect(id_Tcells,id_CD69cells))) %>%
  write_tsv("../out_large/scTenifoldNet/table/IMMUNE_CCA_martina_Net_trimm_removeTCD69cells.tsv")
