# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","Human_Gene_Atlas","Azimuth_Cell_Types_2021")

# GENE SELECTION ----------------------------------------------------------
# list of genes to consider for the enrichment analysis
test_out <- readRDS("out/object/IMMUNE_CCA_martina_trimm_KO_BTK.rds")
glimpse(test_out)

# read in the metadata for the genes
# meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")

# save the ranked object, also change the genes into genenames
full_df <- test_out$diffRegulation %>%
  # left_join(meta,by = c("gene" = "gene_id")) %>%
  mutate(dataset = "test_IMMUNE") %>%
  mutate(log2FC = log2(FC))

list_res_tot <- list(IMMUNE_C1QB = full_df %>%
                       filter(p.adj < 0.05) %>%
                       pull(gene))

# query -------------------------------------------------------------------
list <- lapply(list_res_tot,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  #
  out_enrich
  })

df_enrichr_annotation_enriched_tot <- list$IMMUNE_C1QB %>%
  bind_rows(.id = "annotation")

df_enrichr_annotation_enriched_tot %>%
  write_tsv("out/table/enrichR_out_IMMUNE_BTK.tsv")

library(scales)
df_enrichr_annotation_enriched_tot %>%
  group_by(annotation) %>%
  arrange(P.value) %>%
  dplyr::slice(1:20) %>%
  mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
  mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
  ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,ncol = 1,scales = "free_y")+theme_bw() +
  scale_color_gradientn(colors = c("red","blue"),
                        values = rescale(c(0,1)),
                        limits = c(0,0.2))
  # scale_color_gradient(low = "red",high = "blue")

ggsave("out/image/enrichR_out_IMMUNE_BTK_plot.pdf",width = 7,height = 15)


# shorter versino for martina ---------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020")

# GENE SELECTION ----------------------------------------------------------
# list of genes to consider for the enrichment analysis
test_out <- readRDS("out/object/IMMUNE_CCA_martina_trimm_KO_BTK.rds")
glimpse(test_out)

# read in the metadata for the genes
# meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")

# save the ranked object, also change the genes into genenames
full_df <- test_out$diffRegulation %>%
  # left_join(meta,by = c("gene" = "gene_id")) %>%
  mutate(dataset = "test_IMMUNE") %>%
  mutate(log2FC = log2(FC))

list_res_tot <- list(IMMUNE_C1QB = full_df %>%
                       filter(p.adj < 0.05) %>%
                       pull(gene))

# query -------------------------------------------------------------------
list <- lapply(list_res_tot,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  #
  out_enrich
})

df_enrichr_annotation_enriched_tot <- list$IMMUNE_C1QB %>%
  bind_rows(.id = "annotation")

library(scales)
df_enrichr_annotation_enriched_tot %>%
  group_by(annotation) %>%
  arrange(P.value) %>%
  dplyr::slice(1:20) %>%
  mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
  mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
  ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,ncol = 1,scales = "free")+theme_bw() +
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) + 
  scale_color_gradientn(colors = c("red","blue"),
                        values = rescale(c(0,1)),
                        limits = c(0,0.2)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_discrete(expand = expansion(mult = 0.1))
# scale_color_gradient(low = "red",high = "blue")

ggsave("out/image/enrichR_out_IMMUNE_BTK_plot_short.pdf",width = 7,height = 7)


