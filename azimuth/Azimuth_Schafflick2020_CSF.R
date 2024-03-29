# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(tidyverse)
library(ComplexHeatmap)

# function definition -----------------------------------------------------
#' Transform an NN index
#'
#' @param object Seurat object
#' @param meta.data Metadata
#' @param neighbor.slot Name of Neighbor slot
#' @param key Column of metadata to use
#'
#' @return \code{object} with transfomed neighbor.slot
#'
#' @importFrom SeuratObject Indices
#'
#' @keywords internal
#'
NNTransform <- function(
    object,
    meta.data,
    neighbor.slot = "query_ref.nn",
    key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

# load the data -----------------------------------------------------------
# Download the Azimuth reference and extract the archive
reference <- readRDS(file = "data/human_PBMC_v1.0.0.rds")

# # Load the query object for mapping
object <- readRDS("data/schafflick_csf_harmony_by_patho.rds")
DimPlot(object,reduction = "umap",label = T)
ggsave("out/image/UMAP_schafflickCsf_tot.pdf",width = 7,height = 5)

query <- LoadFileInput(path = "out/object/diet_schafflick_csf.rds")
head(query@meta.data)

# preprocessing -----------------------------------------------------------
# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# label transfer ----------------------------------------------------------
refdata <- lapply(X = c("celltype.l1", "celltype.l2", "celltype.l3"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("celltype.l1", "celltype.l2", "celltype.l3")

query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)

# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

# save the metadata for the query to save the label transfer
query@meta.data %>%
  rownames_to_column("barcodes") %>%
  write_tsv("out/table/schafflick_csf_metaAzimuth.tsv")

query@reductions$proj.umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcodes") %>%
  write_tsv("out/table/schafflick_csf_CoordUMAP_Azimuth.tsv")

# save the query as an object 
saveRDS(query ,file = "out/object/schafflick_csf_queryAzimuth.rds")

# visualization -----------------------------------------------------------
# First predicted metadata field, change to visualize other predicted metadata
id <- c("celltype.l1", "celltype.l2", "celltype.l3")[2]
predicted.id <- paste0("predicted.", id)

# DimPlot of the reference
DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()

# costume visualization
# load the original UMAP in this case
UMAP_query <- read_tsv("out/table/schafflick_csf_CoordUMAP_curated.tsv")
UMAP_query %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point()

meta_query <- read_tsv("out/table/schafflick_csf_metaAzimuth.tsv") %>%
  mutate(robust_score_celltype.l1 = case_when(predicted.celltype.l1.score>0.70&mapping.score>0.70~predicted.celltype.l1,
                                              T~"uncertain"),
         robust_score_celltype.l2 = case_when(predicted.celltype.l2.score>0.70&mapping.score>0.70~predicted.celltype.l2,
                                              T~"uncertain"),
         robust_score_celltype.l3 = case_when(predicted.celltype.l3.score>0.70&mapping.score>0.70~predicted.celltype.l3,
                                              T~"uncertain"))

# define the levels of the cluster variable
level_celltype.l1 <- meta_query %>%
  group_by(predicted.celltype.l1) %>%
  summarise(med = median(predicted.celltype.l1.score)) %>%
  mutate(predicted.celltype.l1 = fct_reorder(predicted.celltype.l1, med,.desc = T)) %>%
  pull(predicted.celltype.l1) %>%
  levels()

level_celltype.l2 <- meta_query %>%
  group_by(predicted.celltype.l2) %>%
  summarise(med = median(predicted.celltype.l2.score)) %>%
  mutate(predicted.celltype.l2 = fct_reorder(predicted.celltype.l2, med,.desc = T)) %>%
  pull(predicted.celltype.l2) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.celltype.l1 = factor(predicted.celltype.l1, levels = level_celltype.l1)) %>%
  ggplot(aes(x=predicted.celltype.l1,y=predicted.celltype.l1.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.70,col="red")
ggsave("out/image/schafflick_csf_predicted.celltype.l1_070.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, levels = level_celltype.l2)) %>%
  ggplot(aes(x=predicted.celltype.l2,y=predicted.celltype.l2.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.70,col="red")
ggsave("out/image/schafflick_csf_predicted.celltype.l2_070.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.celltype.l1 = factor(predicted.celltype.l1, levels = level_celltype.l1)) %>%
  ggplot(aes(y=predicted.celltype.l1,x=predicted.celltype.l1.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.70,col="red",linetype="dashed")
ggsave("out/image/schafflick_csf_predicted.celltype.l1_ridges_070.pdf",height = 5,width = 4)

meta_query %>%
  mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, levels = level_celltype.l2)) %>%
  ggplot(aes(y=predicted.celltype.l2,x=predicted.celltype.l2.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.70,col="red",linetype="dashed")
ggsave("out/image/schafflick_csf_predicted.celltype.l2_ridges_070.pdf",height = 5,width = 4)

# identifyt he most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores
prop_predicted.celltype.l1 <- meta_query %>%
  group_by(seurat_clusters,predicted.celltype.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l1")

prop_predicted2.celltype.l1 <- meta_query %>%
  group_by(predicted.celltype.l1,seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l1")

prop_predicted.celltype.l2 <- meta_query %>%
  group_by(seurat_clusters,predicted.celltype.l2) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l2")

prop_predicted2.celltype.l2 <- meta_query %>%
  group_by(predicted.celltype.l2,seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l2")

pdf("out/image/schafflick_csf_predicted.celltype.l1_heatmapAll.pdf",height = 3,width = 5)
Heatmap(prop_predicted.celltype.l1,
        name = "prop",
        column_title = "celltype l1 all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l1_heatmapAll2.pdf",height = 3,width = 5)
Heatmap(prop_predicted2.celltype.l1,
        name = "prop",
        column_title = "celltype l1 all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l2_heatmapAll.pdf",height = 6,width = 5)
Heatmap(prop_predicted.celltype.l2,
        name = "prop",
        column_title = "celltype l2 all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l2_heatmapAll2.pdf",height = 6,width = 5)
Heatmap(prop_predicted2.celltype.l2,
        name = "prop",
        column_title = "celltype l2 all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# first using all the subcluster annotation, filtering for threshold of scores
prop_predicted.celltype.l1_filter <- meta_query %>%
  filter(robust_score_celltype.l1 != "uncertain") %>%
  group_by(seurat_clusters,predicted.celltype.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l1") %>% 
  as.matrix()

prop_predicted2.celltype.l1_filter <- meta_query %>%
  filter(robust_score_celltype.l1 != "uncertain") %>%
  group_by(predicted.celltype.l1,seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l1") %>% 
  as.matrix()

prop_predicted.celltype.l2_filter <- meta_query %>%
  filter(robust_score_celltype.l2 != "uncertain") %>%
  group_by(seurat_clusters,predicted.celltype.l2) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l2") %>% 
  as.matrix()

prop_predicted2.celltype.l2_filter <- meta_query %>%
  filter(robust_score_celltype.l2 != "uncertain") %>%
  group_by(predicted.celltype.l2,seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l2") %>% 
  as.matrix()

pdf("out/image/schafflick_csf_predicted.celltype.l1_heatmapFilter_070.pdf",height = 3,width = 5)
Heatmap(prop_predicted.celltype.l1_filter,
        name = "prop",
        column_title = "celltype l1 high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l1_heatmapFilter2_070.pdf",height = 3,width = 5)
Heatmap(prop_predicted2.celltype.l1_filter,
        name = "prop",
        column_title = "celltype l1 high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l2_heatmapFilter_070.pdf",height = 6,width = 5)
Heatmap(prop_predicted.celltype.l2_filter,
        name = "prop",
        column_title = "celltype l2 high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("out/image/schafflick_csf_predicted.celltype.l2_heatmapFilter2_070.pdf",height = 6,width = 5)
Heatmap(prop_predicted2.celltype.l2_filter,
        name = "prop",
        column_title = "celltype l2 high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,"barcodes")

# average the position of the clusters
data_query_avg.l1_all <- data_query %>% group_by(predicted.celltype.l1) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

data_query_avg.l2_all <- data_query %>% group_by(predicted.celltype.l2) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot withnout confidence information
ggplot(label= TRUE)+
  geom_point(data = data_query,aes(x = UMAP_1,y = UMAP_2, col = predicted.celltype.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg.l1_all,aes(x = UMAP_1,y = UMAP_2,label = predicted.celltype.l1),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_predicted.celltype.l1_UMAP_all.pdf",width = 7,height = 5)

ggplot(label= TRUE)+
  geom_point(data = data_query,aes(x = UMAP_1,y = UMAP_2, col = predicted.celltype.l2),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg.l2_all,aes(x = UMAP_1,y = UMAP_2,label = predicted.celltype.l2),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_predicted.celltype.l2_UMAP_all.pdf",width = 8,height = 5)

# divide the dataset into uncertain and not
data_query_unc.l1 <- data_query %>%
  filter(robust_score_celltype.l1 == "uncertain")
data_query_unc.l2 <- data_query %>%
  filter(robust_score_celltype.l2 == "uncertain")
#
data_query_defined.l1 <- data_query %>%
  filter(robust_score_celltype.l1 != "uncertain")
data_query_defined.l2 <- data_query %>%
  filter(robust_score_celltype.l2 != "uncertain")

# average the position of the clusters
data_query_avg.l1 <- data_query_defined.l1 %>% group_by(robust_score_celltype.l1) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

data_query_avg.l2 <- data_query_defined.l2 %>% group_by(robust_score_celltype.l2) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc.l1,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined.l1,aes(x = UMAP_1,y = UMAP_2, col = robust_score_celltype.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg.l1,aes(x = UMAP_1,y = UMAP_2,label = robust_score_celltype.l1),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_predicted.celltype.l1_UMAP_robust_070.pdf",width = 7,height = 5)

ggplot(label= TRUE)+
  geom_point(data = data_query_unc.l2,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined.l2,aes(x = UMAP_1,y = UMAP_2, col = robust_score_celltype.l2),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg.l2,aes(x = UMAP_1,y = UMAP_2,label = robust_score_celltype.l2),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_predicted.celltype.l2_UMAP_robust_070.pdf",width = 8,height = 5)

# costume visualization
UMAP_ref <- read_tsv("data/human_PBMC_v1.0.0_UMAPCoorReferenceAzimuth.tsv")
# load the new umap coordiantes from the integration of the datasets
UMAP_ref_query <- read_tsv("out/table/schafflick_csf_CoordUMAP_Azimuth.tsv")
meta_ref <- read_tsv("data/human_PBMC_v1.0.0_metaAzimuth.tsv")
meta_ref_query <- read_tsv("out/table/schafflick_csf_metaAzimuth.tsv") %>%
  mutate(robust_score_celltype.l1 = case_when(predicted.celltype.l1.score>0.70&mapping.score>0.70~predicted.celltype.l1,
                                              T~"uncertain"),
         robust_score_celltype.l2 = case_when(predicted.celltype.l2.score>0.70&mapping.score>0.70~predicted.celltype.l2,
                                              T~"uncertain"),
         robust_score_celltype.l3 = case_when(predicted.celltype.l3.score>0.70&mapping.score>0.70~predicted.celltype.l3,
                                              T~"uncertain"))

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcodes")
data_ref_query <- left_join(UMAP_ref_query,meta_ref_query,"barcodes") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg.l1 <- data_ref %>% group_by(celltype.l1) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)
data_ref_avg.l2 <- data_ref %>% group_by(celltype.l2) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# plor the reference dataset
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2,col=celltype.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg.l1,aes(x = refumap_1,y = refumap_2,label = celltype.l1),col="black")+theme_bw()
# ggsave("out/image/refPBMC_predicted.celltype.l1_UMAP.pdf",width = 7,height = 5)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2,col=celltype.l2),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol = 2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg.l2,aes(x = refumap_1,y = refumap_2,label = celltype.l2),col="black")+theme_bw()
# ggsave("out/image/refPBMC_predicted.celltype.l2_UMAP.pdf",width = 8,height = 5)

# average the position of the clusters
data_ref_query_avg.l1 <- data_ref_query %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref_query, aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_query_avg.l1,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_refPBMC_immune_predicted.celltype.l1_UMAP_all.pdf",width = 7,height = 5)

# average the position of the clusters
data_ref_query_avg2.l1 <- data_ref_query %>%
  filter(robust_score_celltype.l1 != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref_query %>%
               filter(robust_score_celltype.l1 != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_query_avg2.l1,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()
ggsave("out/image/schafflick_csf_refPBMC_immune_predicted.celltype.l1_UMAP_robust_070.pdf",width = 7,height = 5)
