# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
test_microglia <- readRDS("../data/all20_immune.rds")

# wrangling ---------------------------------------------------------------
# extract the metadata
meta <- test_microglia@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(barcode = str_replace(barcode,pattern = "-",replacement = "\\."))

# save the metadata for the astro
write_tsv(meta,file = "data/IMMUNE_martina_meta.tsv")

# define the barcode to filter to
id_barcodes <- meta %>%
  filter(pathology == "c_chronic_active")

# summarise the proportions  
id_barcodes %>%
  group_by(seurat_clusters,disease) %>%
  summarise(n = n())
#
dim(id_barcodes)

# extract the table of counts
test_microglia_raw_data <- test_microglia@assays$RNA@counts %>%
  data.frame() %>%
  dplyr::select(id_barcodes$barcode) %>%
  rownames_to_column(var = "gene")
#
test_microglia_norm_data <- test_microglia@assays$RNA@data %>%
  data.frame() %>%
  dplyr::select(id_barcodes$barcode) %>%
  rownames_to_column(var = "gene")
#
dim(meta)
dim(test_microglia_raw_data)
dim(test_microglia_norm_data)
test_microglia_raw_data[1:10,1:10]
test_microglia_norm_data[1:10,1:10]

# try to trimm some genes to reduce the computation
# distribution of genes that have non negative cells
data.frame(x = rowSums(test_microglia_raw_data[,-1]>0)) %>%
  ggplot(aes(x=x))+geom_histogram()+scale_x_log10()+theme_bw()
#
data.frame(x = rowSums(test_microglia_norm_data[,-1]>0)) %>%
  ggplot(aes(x=x))+geom_histogram()+scale_x_log10()+theme_bw()

# based on the number of cells per clusters I decided to remove the genes that are expressed in less than 18 cells (cluster 10 has 19 cells).
# the idea is to reduce the number of rows that are not bringing information to the generation of the GCN
# remaning genes
sum(rowSums(test_microglia_raw_data[,-1]>0) > 18)
sum(rowSums(test_microglia_norm_data[,-1]>0) > 18)

id_trimming <- rowSums(test_microglia_raw_data[,-1]>0) > 18
id_trimming_norm <- rowSums(test_microglia_norm_data[,-1]>0) > 18

# save the table ----------------------------------------------------------
# all the cells
dim(test_microglia_raw_data)
test_microglia_raw_data %>%
  write_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_raw_data_all.tsv")
#
test_microglia_norm_data %>%
  write_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_norm_data_all.tsv")

test_microglia_raw_data %>%
  dplyr::filter(id_trimming) %>%
  dim()

test_microglia_norm_data %>%
  dplyr::filter(id_trimming_norm) %>%
  dim()

test_microglia_raw_data %>%
  dplyr::filter(id_trimming) %>%
  write_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_raw_data_all_trimm.tsv")

test_microglia_norm_data %>%
  dplyr::filter(id_trimming_norm) %>%
  write_tsv("../out_large/scTenifoldKnk/table/IMMUNE_CCA_martina_norm_data_all_trimm.tsv")
