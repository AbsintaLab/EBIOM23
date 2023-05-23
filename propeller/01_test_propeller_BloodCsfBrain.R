# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)

# read in the data --------------------------------------------------------
# brain
df_tissue_fix <- read_tsv(file = "out/table/meta_azimuth_dfTissue.tsv") %>% 
  # add a broader annotataion for the disease state
  mutate(disease = case_when(pathology_fix=="control"~"brain_control",
                             T~"brain_MS"))

# blood and csf
df_blood_fix <- read_tsv(file = "out/table/meta_azimuth_dfBloodCsf.tsv")

# explore the data --------------------------------------------------------
# confirm the numbers from the tissue dataset
df_tissue_fix %>% 
  group_by(robust_score_celltype.l2) %>% 
  summarise(n = n())

df_blood_fix %>% 
  group_by(robust_score_celltype.l2) %>% 
  summarise(n = n())

# using the predicted score
df_tissue_fix %>% 
  group_by(predicted.celltype.l2) %>% 
  summarise(n = n()) %>% 
  print(n=30)

df_blood_fix %>% 
  group_by(predicted.celltype.l2) %>% 
  summarise(n = n()) %>% 
  print(n=30)

# group by sample
df_tissue_fix %>% 
  filter(robust_score_celltype.l2!="uncertain") %>% 
  group_by(pathology_fix) %>% 
  summarise(n = n())

df_tissue_fix %>% 
  filter(robust_score_celltype.l2!="uncertain") %>% 
  group_by(disease) %>% 
  summarise(n = n())

df_blood_fix %>% 
  filter(robust_score_celltype.l2!="uncertain") %>% 
  group_by(pathology) %>% 
  summarise(n = n())

# using the predicted score not robust
df_tissue_fix %>% 
  filter(predicted.celltype.l2!="uncertain") %>% 
  group_by(pathology_fix) %>% 
  summarise(n = n())

df_tissue_fix %>% 
  filter(predicted.celltype.l2!="uncertain") %>% 
  group_by(disease) %>% 
  summarise(n = n())

df_blood_fix %>% 
  filter(predicted.celltype.l2!="uncertain") %>% 
  group_by(pathology) %>% 
  summarise(n = n())

# run propeller on our metadata -------------------------------------------

# brain MS vs control -----------------------------------------------------
# run it on the predictred l2, not the robust one to avoid the missing values
# Run propeller testing for cell type proportion differences between the groups
df_tissue_fix_filter <- df_tissue_fix %>% 
  filter(!(predicted.celltype.l2 %in% c("Platelet","Eryth","HSPC","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono")))

table(df_tissue_fix_filter$sample_fix,
      df_tissue_fix_filter$dataset)

table(df_tissue_fix_filter$predicted.celltype.l2,
      df_tissue_fix_filter$sample_fix)

out_brain <- propeller(clusters = df_tissue_fix_filter$predicted.celltype.l2,
                       sample = df_tissue_fix_filter$sample_fix,
                       group = df_tissue_fix_filter$disease)

out_brain %>% 
  rownames_to_column("cluster") %>% 
  write_tsv("out/table/propeller_out_tissue_brainMSvsControl.tsv")

# blood MS vs CTRL --------------------------------------------------------
# run it on the predictred l2, not the robust one to avoid the missing values
# Run propeller testing for cell type proportion differences between the groups
df_blood_fix_filter <- df_blood_fix %>% 
  filter(str_detect(pathology,pattern = "^blood")) %>% 
  filter(!(predicted.celltype.l2 %in% c("Platelet","Eryth","HSPC","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono")))

table(df_blood_fix_filter$sample_fix,
      df_blood_fix_filter$pathology)

table(df_blood_fix_filter$predicted.celltype.l2,
      df_blood_fix_filter$sample_fix)

out_blood <- propeller(clusters = df_blood_fix_filter$predicted.celltype.l2,
                       sample = df_blood_fix_filter$sample_fix,
                       group = df_blood_fix_filter$pathology)

out_blood %>% 
  rownames_to_column("cluster") %>% 
  write_tsv("out/table/propeller_out_tissue_bloodMSvsControl.tsv")

# csf MS vs CTRL ----------------------------------------------------------
# run it on the predictred l2, not the robust one to avoid the missing values
# Run propeller testing for cell type proportion differences between the groups
df_CSF_fix_filter <- df_blood_fix %>% 
  filter(str_detect(pathology,pattern = "^CSF")) %>% 
  filter(!(predicted.celltype.l2 %in% c("Platelet","Eryth","HSPC","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono")))

table(df_CSF_fix_filter$sample_fix,
      df_CSF_fix_filter$pathology)

table(df_CSF_fix_filter$predicted.celltype.l2,
      df_CSF_fix_filter$sample_fix)

out_CSF <- propeller(clusters = df_CSF_fix_filter$predicted.celltype.l2,
                     sample = df_CSF_fix_filter$sample_fix,
                     group = df_CSF_fix_filter$pathology)

out_CSF %>% 
  rownames_to_column("cluster") %>% 
  write_tsv("out/table/propeller_out_tissue_csfMSvsControl.tsv")

# all together ------------------------------------------------------------
df_tot_fix_filter <- bind_rows(list(df_CSF_fix_filter %>%
                                      dplyr::select(predicted.celltype.l2,sample_fix,pathology),
                                    df_blood_fix_filter %>% 
                                      dplyr::select(predicted.celltype.l2,sample_fix,pathology),
                                    df_tissue_fix_filter %>% 
                                      dplyr::select(predicted.celltype.l2,sample_fix,pathology = disease)))

table(df_tot_fix_filter$sample_fix,
      df_tot_fix_filter$pathology)

table(df_tot_fix_filter$predicted.celltype.l2,
      df_tot_fix_filter$sample_fix)

out_tot <- propeller(clusters = df_tot_fix_filter$predicted.celltype.l2,
                     sample = df_tot_fix_filter$sample_fix,
                     group = df_tot_fix_filter$pathology)

out_tot %>% 
  rownames_to_column("cluster") %>% 
  write_tsv("out/table/propeller_out_tissue_totMSvsControl.tsv")
