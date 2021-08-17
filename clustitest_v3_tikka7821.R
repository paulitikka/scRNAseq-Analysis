#ClustifyR program test, 15.7.2021
#https://rnabioco.github.io/clustifyr/articles/clustifyR.html
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(sleepwalk)
library(MASS)
library(gplots)
library(SCINA)
library(hdf5r)
library(stringi)
library(limma)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(Biobase)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(Matrix)
library(stats4)
library(BiocGenerics)
library(parallel)
library(DESeq2) #requirese many loadings..
library(usethis)
library(devtools)
library(reshape2)
library(PRROC)
library(WriteXLS)
library(rpart)
library(rpart.plot)
# library(IKAP)
library(stringr)
library(rlist)
require(rlist)
library(gdata)
library(splines)
library(VGAM)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library("xlsx")
library(tidyr)
library(clustifyr) #this is the main program in this file
library(ggplot2)
library(cowplot)
pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1tsne_1721.rds") #you need later..
pbmc.big1 <- FindClusters(pbmc.big1, resolution = 0.55)
pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2tsne_1721.rds") #you need later..
pbmc.big2 <- FindClusters(pbmc.big2, resolution = 0.6) #you need to do this, since the latest is 'seurat_cluster' #This should be 0.7
# Matrix of normalized single-cell RNA-seq counts
# pbmc_matrix <- pbmc.big1$nFeature_RNA
pbmc_matrix=t(FetchData(object = pbmc.big1, vars = rownames(pbmc.big1)))
#default looks like:
# pbmc_matrix <- clustifyr::pbmc_matrix_small
# > pbmc_matrix[0:5,0:5] 
# 5 x 5 sparse Matrix of class "dgCMatrix"
# AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC AAACCGTGCTTCCG AAACCGTGTATGCG
# PPBP         .              .              .              1.566387       .       
# LYZ          1.635873       1.962726       1.995416       4.521175       .       
# S100A9       .              .              .              3.838845       .       
# IGLL5        .              .              .              .              .       
# GNLY         .              .              1.429744       .              3.453545
# > pbmc_matrix[1,1]
# [1] 0
# > pbmc_matrix[1,4] 
# [1] 1.566387
# rownames(pbmc.big1[[]])
# https://stackoverflow.com/questions/42279766/add-metadata-to-seurat-object
# CellsMeta = pbmc.big1@meta.data
# head(CellsMeta)
randomnumbers <- runif(19590, 0.0, 1.1)
pbmc.big1@meta.data["classified"] <- randomnumbers
head(pbmc.big1@meta.data, n=20)
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==0,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==1,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==2,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==3,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==4,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==5,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==6,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==7,]$classified='Granulocytes'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==8,]$classified='Mesenchyme'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==9,]$classified='Endothelial Cells'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==10,]$classified='Dental Epithelia'
pbmc.big1@meta.data[pbmc.big1@meta.data$seurat_clusters==11,]$classified='Mesenchyme'
# pbmc_meta <- clustifyr::pbmc_meta
# https://satijalab.org/seurat/articles/essential_commands.html
vargenes <- VariableFeatures(object = pbmc.big1)
#And the reference..
new_ref_matrix_v3 <- seurat_ref(
  seurat_object = pbmc.big1,        # SeuratV3 object
  cluster_col = "seurat_clusters")    # name of column in meta.data containing cell identities
# metadata: a meta.data table containing the cluster assignments for each cell (not required if a Seurat object is given), otherwise to res as:
# metadata = pbmc_meta, # meta.data table containing cell clusters
# pbmc_meta <- pbmc.big1@meta.data[c('seurat_clusters','classified')]
pbmc.big1@meta.data["tSNE_1"] <- Embeddings(object = pbmc.big1, reduction = "tsne")[,1] 
#if works with tsne, change UMAP_1 to TSNE_1
pbmc.big1@meta.data["tSNE_2"] <- Embeddings(object = pbmc.big1, reduction = "tsne")[,2]
pbmc.big1@meta.data["UMAP_1"] <- Embeddings(object = pbmc.big1, reduction = "tsne")[,1] 
pbmc.big1@meta.data["UMAP_2"] <- Embeddings(object = pbmc.big1, reduction = "tsne")[,2]
pbmc_meta=pbmc.big1@meta.data #do this at the same time.. # Embeddings(object = pbmc.big1, reduction = "tsne")

#For matrix 2, pbmc.big2
# dim(pbmc_matrix2)
pbmc_matrix=c()
pbmc_matrix=t(FetchData(object = pbmc.big2, vars = rownames(pbmc.big2)))
randomnumbers <- runif(18174, 0.0, 1.1)
pbmc.big2@meta.data["classified"] <- randomnumbers
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==0,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==1,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==2,]$classified='Dental Epithelia'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==3,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==4,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==5,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==6,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==7,]$classified='Mesenchyme' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==8,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==9,]$classified='Endothelial Cells' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==10,]$classified='Mesenchyme' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==11,]$classified='Mesenchyme'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==12,]$classified='Ameloblast type 1'
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==13,]$classified='Leukocytes' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==14,]$classified='Stratum intermedium/stellate reticulum'  
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==15,]$classified='Ameloblast type 2' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==16,]$classified='Odontoblast' 
pbmc.big2@meta.data[pbmc.big2@meta.data$seurat_clusters==17,]$classified='Mesenchyme' 
# head(pbmc.big2@meta.data, n=2)
pbmc.big2@meta.data["tSNE_1"] <- Embeddings(object = pbmc.big2, reduction = "tsne")[,1] 
pbmc.big2@meta.data["tSNE_2"] <- Embeddings(object = pbmc.big2, reduction = "tsne")[,2]
pbmc.big2@meta.data["UMAP_1"] <- Embeddings(object = pbmc.big2, reduction = "tsne")[,1] 
pbmc.big2@meta.data["UMAP_2"] <- Embeddings(object = pbmc.big2, reduction = "tsne")[,2]
pbmc_meta=pbmc.big2@meta.data #do this at the same time..
vargenes <- VariableFeatures(object = pbmc.big2)
new_ref_matrix_v3 <- seurat_ref(
  seurat_object = pbmc.big2,        # SeuratV3 object, and name of column in meta.data containing cell identities:
  cluster_col = "seurat_clusters")
plot_dims(
  x = "tSNE_1", # name of column in the meta.data containing the data to plot on x-axis
  y = "tSNE_2", # name of column in the meta.data containing the data to plot on y-axis
  data = pbmc_meta, # meta.data table containing cluster assignments for each cell
  feature = "seurat_clusters") # name of column in meta.data to color cells by,
plot_gene(
  x = "tSNE_1", # name of column in the meta.data containing the data to plot on x-axis
  y = "tSNE_2", # name of column in the meta.data containing the data to plot on y-axis
  expr_mat = pbmc_matrix, # matrix of normalized single-cell RNA-seq counts
  metadata = pbmc_meta %>% rownames_to_column("rn"), # meta.data table containing cluster assignments for each cell
  genes = c("Pitx2", "Fn1"), # vector of gene names to color cells
  cell_col = "rn" )# name of column in meta.data containing the cell IDs

#Here are the calculation procedures..
res=c()
res <- clustify(
  input = pbmc_matrix, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata= pbmc_meta,
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  ref_mat = new_ref_matrix_v3, # matrix of RNA-seq expression data for each cell type, and list of highly variable genes identified with Seurat
  query_genes = vargenes)

res[ order(row.names(res)), ]
# https://stackoverflow.com/questions/32593434/r-change-row-order
res=res[c(1,2,5,6,7,8,9,10,11,12,3,4),]
# res=res[c(1,2,11,12,13,14,15,16,17,18,3,4,5,6,7,8,9,10),] #for the second
plot_cor_heatmap(cor_mat = res) #needs 'classified' column
plot_cor(
  x = "tSNE_1", # name of column in the meta.data containing the data to plot on x-axis
  y = "tSNE_2", # name of column in the meta.data containing the data to plot on y-axis
  cor_mat = res, # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta, # meta.data table containing cluster assignments for each cell
  data_to_plot = colnames(res)[10:11], # name of cell type(s) to plot correlation coefficients, check this e.g. for dental lamina
  cluster_col = "seurat_clusters") # name of column in meta.data containing cell clusters

res2 <- cor_to_call(
  cor_mat = res,                  # matrix correlation coefficients, and name of column in meta.data containing cell clusters:
  cluster_col = "seurat_clusters")
pbmc_meta2 <- call_to_metadata(
  res = res2,                     # data.frame of called cell type for each cluster
  metadata = pbmc_meta,           # original meta.data table containing cell clusters
  cluster_col = "seurat_clusters") # name of column in meta.data containing cell clusters
# https://rnabioco.github.io/clustifyrdata/articles/visualizations.html

corr_umaps <- plot_cor(
  cor_mat = res,                     # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta,              # meta.data table containing UMAP or tSNE data
  data_to_plot = colnames(res)[10:11], # name of cell type(s) to plot correlation coefficients, same as previously, but two at time
  cluster_col = "seurat_clusters")    # name of column in meta.data containing cell clusters
plot_grid(
  # x = "tSNE_1", # name of column in the meta.data containing the data to plot on x-axis
  # y = "tSNE_2", # name of column in the meta.data containing the data to plot on y-axis
  plotlist = corr_umaps,
  rel_widths = c(0.47, 0.53))

clustifyr_types <- plot_best_call(
  cor_mat = res,          # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta,   # meta.data table containing UMAP or tSNE data
  do_label = FALSE,        # should the feature label be shown on each cluster? If you take this way then ok
  do_legend = TRUE,      # should the legend be shown?
  cluster_col = "seurat_clusters")+ggtitle("Seurat cell clusters")+ guides(color = guide_legend(override.aes = list(size=5)))
# https://stackoverflow.com/questions/15059093/ggplot2-adjust-the-symbol-size-in-legends/15059196 
# Compare clustifyr results with known cell identities
known_types <- plot_dims(
  data = pbmc_meta,       # meta.data table containing UMAP or tSNE data
  feature = "classified", # name of column in meta.data to color clusters by
  do_label = FALSE,        # should the feature label be shown on each cluster?
  do_legend = TRUE)+ggtitle("Known cell types")+guides(color = guide_legend(override.aes = list(size=5)))       
# should the legend be shown? TRUE/FALSE
plot_grid(known_types, clustifyr_types)

cbmc_m2=matrix(nrow = 3, ncol = 13)
# :)
# https://stackoverflow.com/questions/18509527/first-letter-to-upper-case/18509816
for (i in 1:3){
  for (j in 1:13){cbmc_m2[i,j]=stringr::str_to_title(cbmc_m[i,j])}}
# https://www.dummies.com/programming/r/how-to-name-matrix-rows-and-columns-in-r/
rownames(cbmc_m2)= c(1,2,3)
colnames(cbmc_m2)= names(cbmc_m)
# #check your 'cbmc_m'; Available metrics include: "hyper", "jaccard", "spearman", "gsea"
list_res=c()
list_res <- clustify_lists(
  input = pbmc_matrix,             # matrix of normalized single-cell RNA-seq counts
  metadata = pbmc_meta,            # meta.data table containing cell clusters
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  marker = cbmc_m2,                # list of known marker genes
  metric = "pct"  )                 # test to use for assigning cell types
# View as heatmap, or plot_best_call
plot_cor_heatmap(
  cor_mat = list_res,              # matrix of correlation coefficients from clustify_lists()
  cluster_rows = FALSE,            # cluster by row?
  cluster_columns = FALSE,         # cluster by column?
  legend_title = "% expressed"  )   # title of heatmap legend

#And then; https://rnabioco.github.io/clustifyrdata/articles/visualizations.html
plot_best_call(
  cor_mat = res, # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta, # meta.data table containing UMAP or tSNE data
  do_label = FALSE, # should the feature label be shown on each cluster?
  do_legend = TRUE, # should the legend be shown?
  cluster_col = "seurat_clusters")+guides(color = guide_legend(override.aes = list(size=5))) +xlab("tSNE_1")+ylab("tSNE_2")
overcluster_test(
  x_col='tSNE_1',
  y_col='tSNE_2',#this finally worked! :)
  expr = pbmc_matrix, # matrix of normalized single-cell RNA-seq counts
  metadata = pbmc_meta, # meta.data table containing UMAP or tSNE data
  ref_mat = new_ref_matrix_v3, # reference matrix containing bulk RNA-seq data for each cell type
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters,# expand cluster number n-fold for overclustering
  n = 2) +guides(color = guide_legend(override.aes = list(size=5)))+xlab("tSNE_1")+ylab("tSNE_2")

comb_ref <- make_comb_ref(
  ref_mat = new_ref_matrix_v3) # reference matrix containing bulk RNA-seq data for each cell type
# # Peek at the new combined reference
comb_ref[1:5, 1:5]
# Run clustifyr() using the combined reference
comb_res <- clustify(
  input = pbmc_matrix, # matrix of normalized single-cell RNA-seq counts
  metadata = pbmc_meta, # meta.data table containing cluster assignments for each cell
  ref_mat = comb_ref, # reference matrix containing bulk RNA-seq data for each cell type
  query_genes = vargenes, # list of highly varible genes identified with Seurat
  cluster_col = "seurat_clusters") # name of column in meta.data containing cell clusters

# Create tSNE and label clusters with the assigned cell types from the combined reference
plot_best_call(
  cor_mat = comb_res, # matrix of correlation coefficients from clustifyr()
  metadata = pbmc_meta, # meta.data table containing UMAP or tSNE data
  do_label = FALSE, # should the feature label be shown on each cluster?
  do_legend = TRUE, # should the legend be shown?
  cluster_col = "seurat_clusters")+guides(color = guide_legend(override.aes = list(size=5)))+xlab("tSNE_1")+ylab("tSNE_2")  
#This works..