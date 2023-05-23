# libraries ---------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(GSEABase)
library(msigdbr)
library(circlize)
library(scales)
library(ComplexHeatmap)
library(ggrepel)

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

# -------------------------------------------------------------------------
obj_WT <- test_out$tensorNetworks$WT
obj_KO <- test_out$tensorNetworks$KO

mat_WT <- as.matrix(obj_WT)
mat_KO <- as.matrix(obj_KO)

# # str(mat1)
# # mat1[1:10,1:10]
# 
# ht_WT <- Heatmap(mat_WT, show_row_names = F,show_column_names = F, 
#                  name = "GCN", 
#                  column_title = "WT", 
#                  #  cluster_rows = F, 
# )
# 
# ht_KO <- Heatmap(mat_KO,show_row_names = F,show_column_names = F,
#                  name = "GCN", 
#                  column_title = "KO",
#                  #  cluster_rows = F, 
# )
# 
# # pdf("images/heatmap_out_vegnette_IMMUNE_WT.pdf",width = 20,height = 20) 
# # draw(ht_WT,heatmap_legend_side = "left",annotation_legend_side = "left")
# # dev.off()
# # 
# # pdf("images/heatmap_out_vegnette_IMMUNE_KO.pdf",width = 15,height = 15) 
# # draw(ht_KO,heatmap_legend_side = "left",annotation_legend_side = "left")
# # dev.off()


# -------------------------------------------------------------------------
# test_out$manifoldAlignment %>%
#   data.frame() %>%
#   rownames_to_column("gene") %>%
#   ggplot(aes(x=NLMA.1,y=NLMA.2,label = gene))+geom_point() + geom_text_repel()

# -------------------------------------------------------------------------
# explore one of the slot in the output
full_df %>%
  filter(p.adj < 0.05) %>%
  dim()

# save a subset list of the significant
full_df %>%
  # filter(p.adj < 0.05) %>%
  write_tsv("out/table/IMMUNE_BTK_KO_all.tsv")

full_df %>%
  filter(p.adj < 0.05) %>%
  write_tsv("out/table/IMMUNE_BTK_KO_significant.tsv")

# full_df %>%
#   ggplot(aes(x=log2FC))+geom_histogram() +
#   theme_bw()
# 
# full_df %>%
#   ggplot(aes(x=p.value))+geom_histogram() +
#   theme_bw()

# simulate a volcano plot
test_plot <- full_df %>%
  mutate(col = case_when(p.adj<0.05 ~1,
         T~0)) %>%
  mutate(col = factor(col))

test_plot %>% 
  ggplot(aes(x=log2FC,y=-log(p.adj))) + 
  # geom_point() 
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FC,y=-log(p.adj),col=col),alpha=0.05)+ 
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FC,y=-log(p.adj),col=col),alpha=0.8)+ 
  # geom_vline(xintercept = c(-0.3,0.3),col="red",linetype="dashed")+ 
  # geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+ 
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+ 
  geom_text_repel( 
    data = test_plot[test_plot$col==1,], 
    aes(label = gene), 
    size = 2, 
    box.padding = unit(0.35, "lines"), 
    point.padding = unit(0.3, "lines")) + 
  theme_bw()
ggsave("out/image/sim_volcano_IMMUNE_BTK_KO.pdf",width = 5,height = 4)


test_out$diffRegulation %>%
  mutate(log2FC = log2(FC)) %>%
  filter(p.adj<0.05) %>%
  dim()

# # -------------------------------------------------------------------------
# # compare the list of gene from my analysis and the paper
# library(UpSetR) 
# 
# # load the significant from the paper
# paper_df <- read_tsv("data/string_protein_annotations_Trem2.tsv")
# 
# # my analysis
# full_df %>%
#   filter(p.adj<0.05) %>%
#   dim()
# 
# sig_edo <- full_df %>%
#   filter(p.adj<0.05) %>%
#   pull(symbol)
# 
# # set up a list with the genes
# list_genes <- list(sig_paper = unique(paper_df$`#node`),
#                    sig_edo = sig_edo)
# 
# upset(fromList(list_genes), order.by = "freq") 
# 
# # which are the genes only in the paper report?
# id <- setdiff(list_genes$sig_paper,list_genes$sig_edo)
# 
# id_gene <- full_df %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.)) %>%
#   filter(symbol %in% id)
# 
# library(ComplexHeatmap)
# test_m <- full_df %>%
#   mutate(rank = nrow(.):1) %>%
#   mutate(rank2 = 1:nrow(.))
# 
# m = matrix(test_m$rank,ncol = 1)
# ha = rowAnnotation(foo = anno_mark(at = id_gene$rank2, labels = id_gene$symbol))
# hm <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha)
# 
# pdf("images/heatmap_out_Trem2_genes.pdf",width = 2,height = 10) 
# draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()
