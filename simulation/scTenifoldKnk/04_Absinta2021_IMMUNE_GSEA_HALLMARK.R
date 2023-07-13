# libraries ---------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(GSEABase)
library(msigdbr)
library(circlize)
library(scales)
library(ggrepel)
library(patchwork)

# read in the ranked object -----------------------------------------------
test_out <- readRDS("out/object/IMMUNE_CCA_martina_trimm_KO_BTK.rds")
glimpse(test_out)

# save the ranked object, also change the genes into genenames
full_df <- test_out$diffRegulation %>%
  # left_join(meta,by = c("gene" = "gene_id")) %>%
  mutate(dataset = "test_IMMUNE") %>%
  mutate(log2FC = log2(FC))

# prepare the dataset with all the annoration needed ---------------------- 
results <- full_df %>%
  split(.$dataset)

# GSEA -------------------------------------------------------------------- 
list_ranks <- lapply(results, function(x){
  
  x <- filter(x,!is.na(gene)) %>%
    filter(!is.infinite(FC)) %>%
    drop_na(FC) %>%
    # average logFC in case of duplicated genenames
    group_by(gene) %>%
    summarise(avg_distance = mean(distance))
  
  ranks <- setNames(x$avg_distance, x$gene)
  ranks
}) 
glimpse(list_ranks)

head((sort(list_ranks$test_IMMUNE,decreasing = T)))

# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
head(h_gene_sets)

pathways <- split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

# RUN GSEA ----------------------------------------------------------------
set.seed(55)
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  # fgsea(pathways, x, minSize=1, maxSize=500, nperm=1000)
  fgsea(pathways, x, minSize=1, maxSize=500,scoreType="pos")  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

head(df_tables_GSEA_all,n=20) 

# sasve the table
write_tsv(df_tables_GSEA_all,file = "out/table/df_tables_GSEA_test_out_IMMUNE_BTK_KO_HALLMARK.tsv")

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 
names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

str(list_collapsedPathways)

# mainPathways 
list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# check the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "out/table/df_tables_GSEA_test_out_IMMUNE_BTK_KO_HALLMARK_non_redundant.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 9)

# test plot to show the main terms in each dataset
# library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
ggsave("out/image/GSEA_test_out_IMMUNE_BTK_KO_non_redundant_HALLMARK.pdf",width = 5,height = 5)

# -------------------------------------------------------------------------
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  geom_text_repel(size = 2,
                  point.padding = 0, # additional padding around each point
                  min.segment.length = 0, # draw all line segments
                  max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                  box.padding = 0.8,
                  segment.alpha = 0.1,max.overlaps = 8)
ggsave("out/image/GSEA_test_out_IMMUNE_BTK_KO_all_HALLMARK.pdf",width = 5,height = 5)

# plot profile ------------------------------------------------------------
# library("patchwork")
plot_pathway <- lapply(list_mainPathways$test_IMMUNE[1:5],function(x){
  # name <- str_sub(x,start = 1,end = -4)
  name <- x
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
  plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
})

wrap_plots(plot_pathway)
ggsave(filename = "out/image/GSEA_plot_profile_IMMUNE_BTK_KO_non_redundant_HALLMARK.pdf",width = 15,height = 10)

# -------------------------------------------------------------------------
id_pathways <- list_tables_GSEA_all$test_IMMUNE %>%
  dplyr::slice(1:9)

plot_pathway <- lapply(id_pathways$pathway,function(x){
  # name <- str_sub(x,start = 1,end = -4)
  name <- x
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
  plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
})

wrap_plots(plot_pathway)
ggsave(filename = "out/image/GSEA_plot_profile_IMMUNE_BTK_KO_all_HALLMARK.pdf",width = 15,height = 15)
