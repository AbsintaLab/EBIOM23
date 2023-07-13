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

# save the ranked object, also change the genes into genenames
full_df <- test_out$diffRegulation %>%
  # left_join(meta,by = c("gene" = "gene_id")) %>%
  mutate(dataset = "test_IMMUNE") %>%
  mutate(log2FC = log2(FC))

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
