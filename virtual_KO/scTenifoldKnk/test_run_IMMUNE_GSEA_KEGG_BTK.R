# libraries ---------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(GSEABase)
library(msigdbr)
library(circlize)
library(scales)

# read in the ranked object -----------------------------------------------
test_out <- readRDS("out/object/IMMUNE_CCA_martina_trimm_KO_BTK.rds")
glimpse(test_out)

# read in the metadata for the genes
# meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")

# save the ranked object, also change the genes into genenames
full_df <- test_out$diffRegulation %>%
  # left_join(meta,by = c("gene" = "gene_id")) %>%
  mutate(dataset = "test_IMMUNE") %>%
  mutate(log2FC = log2(FC))

# prepare the dataset with all the annoration needed ---------------------- 
results <- full_df %>%
  split(.$dataset)

# # give the names 
# names(results) <- str_remove_all(file,pattern = ".tsv")
# results

# GSEA -------------------------------------------------------------------- 
# use the FC dataset to create the ranked list of genes 
# Symbol or Entrez? 
# filter out the Inf rations
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
# important, for this dataset the genenames used have been the homo sapiens gene names. therefore to avoid issues I have to use the homo sapiens annotation for the enrichemnt analysis
# reactome_gene_sets <- msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP:REACTOME")
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
# reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
# head(reactome_gene_sets)
head(kegg_gene_sets)
# head(reactome_gene_sets)

# format in order to be accepted by GSEA
# pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
pathways <- split(x = kegg_gene_sets$gene_symbol, f = kegg_gene_sets$gs_name)
# pathways <- split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
# pathways <- split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)
# pathways <- split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
# head(pathways)

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

# # generate an heatmap for the comparison ofthe nes across datasets and pathways
# mat_dataset <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,pathway,NES) %>%
#   mutate(pathway = str_remove(pathway,pattern = "_review")) %>%
#   pivot_wider(names_from = dataset,values_from = NES) %>%
#   column_to_rownames(var = "pathway")

#
# library(tidyverse)
# library(circlize)
# library(ComplexHeatmap)
# 
# ht <- Heatmap(mat_dataset,name = "NES", column_title = "GSEA old vs young")
# 
# pdf(file = "images/heatmap_NES_targeted.pdf", width = 8, height = 4)
# draw(ht,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 70), "mm"))
# dev.off()

# dim(df_tables_GSEA_all)

# df_tables_GSEA_all %>%
#   filter(str_detect(pathway,pattern = "ALZH"))

head(df_tables_GSEA_all,n=20) 

# sasve the table
write_tsv(df_tables_GSEA_all,file = "out/table/df_tables_GSEA_test_out_IMMUNE_BTK_KO_KEGG.tsv")
# write_tsv(df_tables_GSEA_all,file = "out/df_tables_GSEA_test_out_IMMUNE_BTK_KO_HALLMARK.tsv")
# write_tsv(df_tables_GSEA_all,file = "out/df_tables_GSEA_test_out_IMMUNE_BTK_KO_REACTOME.tsv")

# plot the rank as an heatmap
# df_tables_GSEA_all %>%
#   mutate(rank = rank(NES)) %>%
#   # filter(NES>0) %>%
#   ggplot(aes(x=dataset,y=rank,fill = NES))+geom_tile()+scale_fill_gradient(low = "blue",high = "red")

# df_tables_GSEA_all %>%
#   mutate(rank = 1:nrow(.)) %>%
#   # mutate(rank = rank(NES)) %>%
#   # filter(NES>0) %>%
#   ggplot(aes(x=dataset,y=rank,fill = rank))+geom_tile() + 
#   scale_fill_gradient(low = "blue",high = "red")
#   # scale_fill_gradientn(colors=c("red","white","blue"))
#   # scale_fill_continuous()
#   # scale_fill_gradientn(colors=c("red","white","blue"),
#   #                      values = rescale(c(2,0,-1)),
#   #                      limits=c(-1,2))

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

# id_pathways <- list_tables_GSEA_all$test_IMMUNE %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.)) %>%
#   filter(str_detect(pathway,pattern = "ALZHEIMER|LYSOSOME|CHOLESTER|OXIDATIVE"))
# 
# library(ComplexHeatmap)
# test_m <- df_tables_GSEA_all %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.))
# 
# m = matrix(test_m$rank,ncol = 1)
# ha = rowAnnotation(foo = anno_mark(at = id_pathways$rank2, labels = id_pathways$pathway))
# hm <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha)
# 
# # pdf("images/heatmap_out_Trem2_NES.pdf",width = 4.5,height = 5) 
# draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()

# explore the correlation between each subset focussing on the ES
# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# # split the dataset per type
# list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# collapsedPathways_MG <- collapsePathways(df_tables_GSEA_all %>%
#                                            filter(dataset %in% "list_df_comparisons_MG_old_young"), pathways, list_ranks$list_df_comparisons_MG_old_young)
# glimpse(collapsedPathways_MG)
str(list_collapsedPathways)

# filter only the pathways in the fgsea object from the collapsed one 
# mainPathways_MG <- df_tables_GSEA_all %>%
#   filter(pathway %in% collapsedPathways$mainPathways) %>%
#   arrange(padj,-abs(NES)) %>%
#   pull(pathway) 
# 
# mainPathways 

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# produce the plot 
# plotGseaTable(pathways[list_mainPathways$res_sen_final_shr],
#               list_ranks$res_sen_final_shr,
#               list_tables_GSEA_all$res_sen_final_shr, gseaParam = 0.5)

# save list of non redundant terms

# chackt the order of the names is the same
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
  write_tsv(file = "out/table/df_tables_GSEA_test_out_IMMUNE_BTK_KEGG_KO_non_redundant.tsv")
# df_tables_GSEA_all_non_redundant %>%
#   write_tsv(file = "out/df_tables_GSEA_test_out_IMMUNE_BTK_KO_HALLMARK_non_redundant.tsv")
# df_tables_GSEA_all_non_redundant %>%
#   write_tsv(file = "out/df_tables_GSEA_test_out_IMMUNE_BTK_KO_REACTOME_non_redundant.tsv")
# df_tables_GSEA_all %>%
#   filter(dataset == "ENDO") %>%
#   filter(pathway %in% list_mainPathways$ENDO) %>%
#   write_tsv(file = "out/table_GSEA_non_redundant_REACTOME.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 9)

# test plot to show the main terms in each dataset
library(ggrepel)
# df_tables_GSEA_all_non_redundant %>%
#   # shorten the label of the pathway
#   mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
#            str_sub(start = 1,end = 35)) %>%
#   # mutate(min_log10_padj = -log10(padj)) %>%
#   ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
#   geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
# ggsave("images/GSEA_test_out_IMMUNE_BTK_KO_non_redundant_REACTOME.pdf",width = 5,height = 5)

df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
ggsave("out/image/GSEA_test_out_IMMUNE_BTK_KO_non_redundant_KEGG.pdf",width = 5,height = 5)

# df_tables_GSEA_all_non_redundant %>%
#   # shorten the label of the pathway
#   mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_") %>%
#            str_sub(start = 1,end = 35)) %>%
#   # mutate(min_log10_padj = -log10(padj)) %>%
#   ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
#   geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
# ggsave("images/GSEA_test_out_IMMUNE_BTK_KO_non_redundant_HALLMARK.pdf",width = 5,height = 5)

# -------------------------------------------------------------------------
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
  geom_text_repel(size = 2,
                  point.padding = 0, # additional padding around each point
                  min.segment.length = 0, # draw all line segments
                  max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                  box.padding = 0.8,
                  segment.alpha = 0.1,max.overlaps = 8)
ggsave("out/image/GSEA_test_out_IMMUNE_BTK_KO_all_KEGG.pdf",width = 5,height = 5)

# df_tables_GSEA_all %>%
#   # shorten the label of the pathway
#   mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_") %>%
#            str_sub(start = 1,end = 35)) %>%
#   # mutate(min_log10_padj = -log10(padj)) %>%
#   ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
#   geom_text_repel(size = 2,
#                   point.padding = 0, # additional padding around each point
#                   min.segment.length = 0, # draw all line segments
#                   max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
#                   box.padding = 0.8,
#                   segment.alpha = 0.1,max.overlaps = 8)
# ggsave("images/GSEA_test_out_IMMUNE_BTK_KO_all_HALLMARK.pdf",width = 5,height = 5)
# 
# df_tables_GSEA_all %>%
#   # shorten the label of the pathway
#   mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
#            str_sub(start = 1,end = 35)) %>%
#   # mutate(min_log10_padj = -log10(padj)) %>%
#   ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() + 
#   geom_text_repel(size = 2,
#                   point.padding = 0, # additional padding around each point
#                   min.segment.length = 0, # draw all line segments
#                   max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
#                   box.padding = 0.8,
#                   segment.alpha = 0.1,max.overlaps = 8)
# ggsave("images/GSEA_test_out_IMMUNE_BTK_KO_all_REACTOME.pdf",width = 5,height = 5)

# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# show an heatmap for hte rank and highlight some of the pathways

# # PLOT PROFILE ------------------------------------------------------------ 
library("patchwork")

plot_pathway <- lapply(list_mainPathways$test_IMMUNE[1:7],function(x){
  # name <- str_sub(x,start = 1,end = -4)
  name <- x
  # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
  plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
})

wrap_plots(plot_pathway)
ggsave(filename = "out/image/GSEA_plot_profile_IMMUNE_BTK_KO_non_redundant_KEGG.pdf",width = 15,height = 10)

#
# plot_pathway <- lapply(list_mainPathways$test_IMMUNE[1:5],function(x){
#   # name <- str_sub(x,start = 1,end = -4)
#   name <- x
#   # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
#   plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
# })
# 
# wrap_plots(plot_pathway)
# ggsave(filename = "images/GSEA_plot_profile_IMMUNE_BTK_KO_non_redundant_HALLMARK.pdf",width = 15,height = 10)
# 
# plot_pathway <- lapply(list_mainPathways$test_IMMUNE[1:9],function(x){
#   # name <- str_sub(x,start = 1,end = -4)
#   name <- x
#   # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
#   plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
# })
# 
# wrap_plots(plot_pathway)
# ggsave(filename = "images/GSEA_plot_profile_IMMUNE_BTK_KO_non_redundant_REACTOME.pdf",width = 15,height = 10)

# select a subset of pathways
# id_pathways <- list_tables_GSEA_all$test_Trem2 %>%
#   filter(str_detect(pathway,pattern = "ALZHEIMER|LYSOSOME|CHOLESTER|OXIDATIVE")) %>%
#   pull(pathway)

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
ggsave(filename = "out/image/GSEA_plot_profile_IMMUNE_BTK_KO_all_KEGG.pdf",width = 15,height = 15)

# plot_pathway <- lapply(id_pathways$pathway,function(x){
#   # name <- str_sub(x,start = 1,end = -4)
#   name <- x
#   # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3)
#   plotEnrichment(pathways[[x]], list_ranks$test_IMMUNE) + labs(title = name)
# })
# 
# wrap_plots(plot_pathway)
# ggsave(filename = "images/GSEA_plot_profile_IMMUNE_BTK_KO_all_HALLMARK.pdf",width = 15,height = 15)
