#Pauli Tikka, scRNAseq Seurat analysis pipeline for 10X data, ver 4.6.2021
#First install Seurat preferably in R studio, you may need to install many things 
#if this is your first or second etc. session.. if you are prompted to install dependencies, please do so! or a in everything
BiocManager::install(c('S4Vectors', 'IRanges', 'GenomicRanges',
                       'GenomeInfoDb', 'Biobase', 'SummarizedExperiment',
                       'MatrixGenerics', 'matrixStats', 'Matrix.utils',
                       'stats4', 'BiocGenerics', 'parallel',
                       'DESeq2', 'scran',
                       'scater', 'xlsx', 'tidyr'))
install.packages("stringi",type="win.binary") #you may need stringi
install.packages('dplyr')
install.packages('Seurat')
install.packages('tidyverse')
install.packages('sleepwalk')
install.packages('SCINA')
install.packages('hdf5r')
install.packages("devtools")
install.packages("SeuratWrappers")
devtools::install_github("NHLBI-BCB/IKAP")
#Installing monocle can be hard if you have R 4.0..
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
remotes::install_github('satijalab/seurat-wrappers') #you need this!
install.packages("namespace")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install(c('Rsamtools'))
BiocManager::install(c('ComplexHeatmap'))
install.packages("Rsamtools")
install.packages("Signac")
# source("IKAP_Seurat3.R")
# getwd()
# setwd("D:/Codes and Instructions")
install.packages("xlsx")
install.packages('Scatter')
BiocManager::install('Scran') #also for ('DESeq2')
install.packages('matrixStats') #for additions
BiocManager::install('matrixStats')
#Droplets away
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)

#1) Then load the packages and data
setwd(dir='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat')
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
library(stats4)
library(parallel)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(Matrix)
library(Biobase)
library(matrixStats)
library(MatrixGenerics)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DESeq2) #requirese many loadings..
library(usethis)
library(devtools)
library(reshape2)
library(PRROC)
library(WriteXLS)
library(rpart)
library(rpart.plot)
set.seed(1234)
# library(IKAP)
library(stringr)
library(rlist)
require(rlist)
library(gdata)
library(splines)
library(VGAM)
library(scuttle)
library(scran)
library(scater)
library("xlsx")
library(tidyr)
library("monocle3")
library(parallel)
library(BiocGenerics)
library(Biobase)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(bit)
library(bit64)
library(matrixStats)
library(stats4)
library(S4Vectors)
library(IRanges)
library(MatrixGenerics)
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(plyr)
library(lattice)
library(Rmisc)
library(monocle3)
library(dplyr) 
library(XVector)
library(Biostrings)
library(Rsamtools)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(shiny)
library(Matrix)
library(ggplot2)
library(patchwork)
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3) #library(rbind)
library(base)
library(Matrix)
library(namespace)
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(cowplot))
# suppressMessages(require(org.Hs.eg.db))
#library() #for addition

#Standard workflow (below link) was :
# https://satijalab.org/seurat/articles/essential_commands.html
# Load the PBMC dataset in case you start from scratch:
# pbmc1.data <- Read10X(data.dir = "D:/Data (Raw) from Tooth Development/10X-data_2021_E11.5_dental lamina/mandible1_count/outs/raw_feature_bc_matrix")
# pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "E11.5_1", min.cells = 3, min.features = 200)
pbmc1.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E11.5_dental lamina/mandible1_count/outs/filtered_feature_bc_matrix.h5")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "E11.5_1", min.cells = 3, min.features = 200)
pbmc2.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E11.5_dental lamina/mandible2_count/outs/filtered_feature_bc_matrix.h5")
pbmc2 <- CreateSeuratObject(counts = pbmc2.data, project = "E11.5_2", min.cells = 3, min.features = 200)
pbmc3.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E11.5_dental lamina/mandible3_count/outs/filtered_feature_bc_matrix.h5")
pbmc3 <- CreateSeuratObject(counts = pbmc3.data, project = "E11.5_3", min.cells = 3, min.features = 200)
pbmc4.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E11.5_dental lamina/mandible4_count/outs/filtered_feature_bc_matrix.h5")
pbmc4 <- CreateSeuratObject(counts = pbmc4.data, project = "E11.5_4", min.cells = 3, min.features = 200)
pbmc1l.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E14.25/1_1_count/outs/filtered_feature_bc_matrix.h5")
pbmc1l <- CreateSeuratObject(counts = pbmc1l.data, project = "E14.25_1", min.cells = 3, min.features = 200)
pbmc2l.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E14.25/1_2_count/outs/filtered_feature_bc_matrix.h5")
pbmc2l <- CreateSeuratObject(counts = pbmc2l.data, project = "E14.25_2", min.cells = 3, min.features = 200)
pbmc3l.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E14.25/1_count/outs/filtered_feature_bc_matrix.h5")
pbmc3l <- CreateSeuratObject(counts = pbmc3l.data, project = "E14.25_3", min.cells = 3, min.features = 200)
pbmc4l.data <- Read10X_h5("D:/Data (Raw) from Tooth Development/10X-data_2021_E14.25/2_count/outs/filtered_feature_bc_matrix.h5")
pbmc4l <- CreateSeuratObject(counts = pbmc4l.data, project = "E14.25_4", min.cells = 3, min.features = 200)
#I rather put everything that leave away
gc()
memory.limit(9999999999)
rm(pbmc1.data,pbmc2.data,pbmc3.data,pbmc4.data,pbmc1l.data,pbmc2l.data,pbmc3l.data,pbmc4l.data)
######
pbmc1$dataset <- 'Initiation1'
pbmc2$dataset <- 'Initiation2'
pbmc3$dataset <- 'Initiation3'
pbmc4$dataset <- 'Initiation4'
pbmc1l$dataset <- 'Later1'
pbmc2l$dataset <- 'Later2'
pbmc3l$dataset <- 'Later3'
pbmc4l$dataset <- 'Later4'
# sce <- CreateSeuratObject(assays = list(counts = cbind(pbmc1, pbmc2, pbmc3, pbmc4)))
#Should you have well processed data:
# pbmc1 <- readRDS(file ="D:/Seurat/filtered_gene_bc_matrices/mandible1_count/pbmc_tooth_tikka26521.rds") #check syntax
# #pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "PBMC_Initiation1", min.cells = 3, min.features = 200) #this is not needed...
# pbmc2 <- readRDS(file ="D:/Seurat/filtered_gene_bc_matrices/mandible2_count/pbmc_tooth_tikka26521.rds")
#######
#Merging Works like this: # https://satijalab.org/signac/articles/merging.html
pbmc.big1 <- merge(pbmc1, y = c(pbmc2, pbmc3, pbmc4), add.cell.ids = c("Initiation1", "Initiation2", "Initiation3", "Initiation4"), project = "Initiation (E11.5)") 
pbmc.big2 <- merge(pbmc1l, y = c(pbmc2l, pbmc3l, pbmc4l), add.cell.ids = c("Later1", "Later2", "Later3", "Later4"), project = "Later (E14.25)")
gc()
memory.limit(9999999999)
rm(pbmc1,pbmc2,pbmc3,pbmc4,pbmc1l,pbmc2l,pbmc3l,pbmc4l)
pbmc.big1[["percent.mt"]] <- PercentageFeatureSet(pbmc.big1, pattern = "^mt-")
pbmc.big2[["percent.mt"]] <- PercentageFeatureSet(pbmc.big2, pattern = "^mt-")
# rm(pbmc.big1,pbmc.big2)
#Alternative way (better) for filtering (if needed...)
apply(
  pbmc.big1@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
) -> pbmc.big1$Percent.Largest.Gene

apply(
  pbmc.big2@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
) -> pbmc.big2$Percent.Largest.Gene # head(pbmc.big1$Percent.Largest.Gene)
#Filtering ribosomal genes:
# https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
rb.genes <- rownames(pbmc.big1)[grep("^Rp[sl]",rownames(pbmc.big1))]
C<-GetAssayData(object = pbmc.big1, slot = "counts")
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
pbmc.big1 <- AddMetaData(pbmc.big1, percent.ribo, col.name = "percent.ribo")
VlnPlot(pbmc.big1, features = "percent.ribo", pt.size = 0.1)+ ggtitle("")+ylim(3,45)+xlab("Reads in Cell")+ylab("Ribosomal Genes in Cell (%)")+ NoLegend()
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html
rb.genes <- rownames(pbmc.big2)[grep("^Rp[sl]",rownames(pbmc.big2))]
C<-GetAssayData(object = pbmc.big2, slot = "counts")
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
pbmc.big2 <- AddMetaData(pbmc.big2, percent.ribo, col.name = "percent.ribo")
VlnPlot(pbmc.big2, features = "percent.ribo", pt.size = 0.1)+ ggtitle("")+ylim(0,4)+xlab("Reads in Cell")+ylab("Ribosomal Genes in Cell (%)")+ NoLegend()

#Let us not take hb (28.6.21, tikka).. let us take them! (1.7.21)
hb.genes <- rownames(pbmc.big1)[grep("^Hb",rownames(pbmc.big1))]
C<-GetAssayData(object = pbmc.big1, slot = "counts")
percent.hb <- colSums(C[hb.genes,])/Matrix::colSums(C)*100
pbmc.big1 <- AddMetaData(pbmc.big1, percent.hb, col.name = "percent.hb")
VlnPlot(pbmc.big1, features = "percent.hb", pt.size = 0.1)+ ggtitle("")+ylim(0,1)+xlab("Reads in Cell")+ylab("Blood Genes in Cell (%)")+ NoLegend()
#0.5%is ok...
hb.genes <- rownames(pbmc.big2)[grep("^Hb",rownames(pbmc.big2))]
C<-GetAssayData(object = pbmc.big2, slot = "counts")
percent.hb <- colSums(C[hb.genes,])/Matrix::colSums(C)*100
pbmc.big2 <- AddMetaData(pbmc.big2, percent.hb, col.name = "percent.hb")
VlnPlot(pbmc.big2, features = "percent.hb", pt.size = 0.1)+ ggtitle("")+ylim(0,0.5)+xlab("Reads in Cell")+ylab("Blood Genes in Cell (%)")+ NoLegend()

#First backup save:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1i_16821.rds")
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1i_16821.rds") #you need later..
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2i_16821.rds")
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2i_16821.rds") #you need later..
#For number of cells per replica or later cluster:
table ( Idents(pbmc.big1) )
table ( Idents(pbmc.big2) )

#From now on, it is better to analyse one experimental type (big1 and big2) separately:
#pbmc.big
#The data should be handeled, but just in case
VlnPlot(pbmc.big1 , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# https://www.biostars.org/p/407036/
# https://github.com/satijalab/seurat/issues/259
plot1 <- FeatureScatter(pbmc.big2 , feature1 = "nCount_RNA", feature2 = "percent.mt") #code changed to match pbmc.big
plot2 <- FeatureScatter(pbmc.big2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ggtitle("")+ylim(0,40)+xlab("Reads in Cell")+ylab("Mitochondrial Genes in Cell (%)")
plot2 + ggtitle("")+xlab("Reads in Cell")+ylab("Genes in Cell")
# rm(plot1, plot2)
#pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 8 ) #Check, or create automatic way to do so. Not constant.
VlnPlot(pbmc.big2,features="nCount_RNA")+ ggtitle("")+xlab("Replica Number")+ylab("Number of Reads")
VlnPlot(pbmc.big2,features="nFeature_RNA")+  ggtitle("")+ylim(300,1000)+xlab("Replica Number")+ylab("Number of Genes")
VlnPlot(pbmc.big2,features="Percent.Largest.Gene")+ ggtitle("")+ylim(0,20)+xlab("Replica Number")+ylab("Percentage of Largest Gene (%)")
VlnPlot(pbmc.big2, features=c("percent.mt"))+ ggtitle("")+ylim(0,10)+xlab("Reads in Cell")+ylab("Mitochondrial Genes in Cell (%)")
#Remove title by adding title that has nothing, i.e. ("")
  #https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Automated_Cell_Type_Annotation
# Iterate max and min number of nFetures with the below equation.
# The equation needs to yield (below) 0.99 and (above) 0.01 with the limit values, e.g. 5500 and 500. 
# For the lower limit (e.g. 500) check the plot for the nFeature value, where percent.mt values 
#  have their last most colorful spots (e.g. 1000), and take a (below) weighted average of the two, e.g.
# (320*3+1800*2)/5. Round up with 100 (or 10 if the value is under 1000):
# length(pbmc.big1$nFeature_RNA[pbmc.big1$nFeature_RNA<5200])/length(pbmc.big1$nFeature_RNA) 

#Initiation experiment (E11.5)
a=quantile(pbmc.big1$nFeature_RNA[ pbmc.big1$nFeature_RNA<quantile(pbmc.big1$nFeature_RNA,0.25) ], p=0.5) #596 #check this limit, it is usually ~1000
b=sd((pbmc.big1$nFeature_RNA[ pbmc.big1$nFeature_RNA<quantile(pbmc.big1$nFeature_RNA,0.25) & !pbmc.big1$nFeature_RNA==0 ])) #179.3488
round(a+b,0) #499->500
a=quantile(pbmc.big1$nFeature_RNA[ pbmc.big1$nFeature_RNA>quantile(pbmc.big1$nFeature_RNA,0.9) ], p=0.5) #596 #check this limit, it is usually ~1000
b=sd((pbmc.big1$nFeature_RNA[ pbmc.big1$nFeature_RNA>quantile(pbmc.big1$nFeature_RNA,0.9) ])) #179.3488
round(a+b,0) #4552 ->4600
VlnPlot(pbmc.big1,features="nFeature_RNA")+  ggtitle("")+ylim(200,1000)+xlab("Replica Number")+ylab("Number of Genes")
VlnPlot(pbmc.big1,features="nFeature_RNA")+  ggtitle("")+xlab("Replica Number")+ylab("Number of Genes")+ylim(1500,6000)
#Later experiment (E14.25)
a=quantile(pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA<quantile(pbmc.big2$nFeature_RNA,0.25) ], p=0.5) #596 #check this limit, it is usually ~1000
b=sd((pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA<quantile(pbmc.big2$nFeature_RNA,0.25) & !pbmc.big2$nFeature_RNA==0 ])) #179.3488
round(a+b,0) #676->680
a=quantile(pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA>quantile(pbmc.big2$nFeature_RNA,0.9) ], p=0.5) #596 #check this limit, it is usually ~1000
b=sd((pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA>quantile(pbmc.big2$nFeature_RNA,0.9) ])) #179.3488
round(a+b,0) #5008->5000
# quantile(pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA>3000 ], p=0.5) #4202 #check this it is usually ~3000
# sd(pbmc.big2$nFeature_RNA[ pbmc.big2$nFeature_RNA>3000 ]) #806.6213
VlnPlot(pbmc.big2,features="nFeature_RNA")+  ggtitle("")+ylim(200,1200)+xlab("Replica Number")+ylab("Number of Genes")
VlnPlot(pbmc.big2,features="nFeature_RNA")+  ggtitle("")+xlab("Replica Number")+ylab("Number of Genes")+ylim(1500,7000)
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html#Plot_QC

#For mitochodrias: Iterate the below x (e.g. start 25) so that you get more than 0.95 or 0.8,(now for pbmc.big1 it was 24)
# length(pbmc.big1$percent.mt[pbmc.big1$percent.mt<15])/length(pbmc.big1$percent.mt)
a=quantile(pbmc.big1$percent.mt, p=0.75)
b=2^sd(log2(pbmc.big1$percent.mt[!pbmc.big1$percent.mt==0]))
round(a+b) #15 ok
a=quantile(pbmc.big2$percent.mt, p=0.75)
b=2^sd(log2(pbmc.big2$percent.mt[!pbmc.big2$percent.mt==0]))
round(a+b) #6 ok
# round(10^quantile(sd(log10(pbmc.big1$percent.mt)))) #7.54 ->5.01+7.54=12.55~12 (round down nearest integer, needs to be stringent)
# quantile(pbmc.big2$percent.mt, p=0.5) #2.22
# sd(pbmc.big2$percent.mt) #3.95 ->2.22+3.95=6.17~6
#For the largest gene:
a=quantile(pbmc.big1$Percent.Largest.Gene, p=0.9)#90% is ok
b=2^sd(log2(pbmc.big1$Percent.Largest.Gene))
round(a+b) #8%
a=quantile(pbmc.big2$Percent.Largest.Gene, p=0.9)#90% is ok
b=2^sd(log2(pbmc.big2$Percent.Largest.Gene))
round(a+b) #11%
# round(10^quantile(log10(pbmc.big1$Percent.Largest.Gene), p=0.5)+10^sd(log10(pbmc.big1$Percent.Largest.Gene)),0) 
#~16 (round normally)
#If skewed distribution the percent.largest.gene should be c, see below, otherwise just 'mean+sd' (and round up)
# a=round(median(log(pbmc.big2$Percent.Largest.Gene, exp(1))),3)
# b=round(sd(log(pbmc.big2$Percent.Largest.Gene, exp(1))),3)
# c=round(exp(a+b),0)+1 
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html#Plot_QC

#For ribos: Iterate the below x (e.g. start 25) so that you get more than 0.95 or 0.8,(now for pbmc.big1 it was 24)
# length(pbmc.big1$percent.mt[pbmc.big1$percent.mt<15])/length(pbmc.big1$percent.mt)
a=quantile(pbmc.big1$percent.ribo, p=0.05)
b=2^sd(log2(pbmc.big1$percent.ribo[!pbmc.big1$percent.ribo==0 & pbmc.big1$percent.ribo<10 ]))
round(a-b) #8, ok also from picture:
VlnPlot(pbmc.big1, features = "percent.ribo", pt.size = 0.1)+ ggtitle("")+ylim(0,40)+xlab("Reads in Cell")+ylab("Ribosomal Genes in Cell (%)")+ NoLegend()
a=quantile(pbmc.big2$percent.ribo, p=0.05)
b=2^sd(log2(pbmc.big2$percent.ribo[!pbmc.big2$percent.ribo==0 & pbmc.big2$percent.ribo<8 ]))
round(a-b) #5 ok, check the image
VlnPlot(pbmc.big2, features = "percent.ribo", pt.size = 0.1)+ ggtitle("")+ylim(0,50)+xlab("Reads in Cell")+ylab("Ribosomal Genes in Cell (%)")+ NoLegend()
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html
#This value can be lower, since now trying to remove the values

#For blood:
a=quantile(pbmc.big1$percent.hb, p=0.75)
# b=2^sd(log2(pbmc.big1$percent.hb[!pbmc.big1$percent.hb==0]))
round(a,2) #0.84, ok also from picture
a=quantile(pbmc.big2$percent.hb, p=0.75)
# b=2^sd(log2(pbmc.big2$percent.ribo[!pbmc.big2$percent.ribo==0]))
round(a,2) #0.23 ok, check the image
VlnPlot(pbmc.big2, features = "percent.hb", pt.size = 0.1)+ ggtitle("")+ylim(0,0.5)+xlab("Reads in Cell")+ylab("Blood Genes in Cell (%)")+ NoLegend()

#For number of cells per replica or later cluster:
table ( Idents(pbmc.big1) )
table ( Idents(pbmc.big2) ) #should look like:
# E11.5_1 E11.5_2 E11.5_3 E11.5_4 
# 7041    7460    8534    8703 
# E14.25_1 E14.25_2 E14.25_3 E14.25_4 
# 3091     5290    11504     7609 
# dim(pbmc.big1)
# pbmc.big1 <- pbmc.big1[!grepl("MALAT1", rownames(pbmc.big1)), ] #MALAT1 is not there..
# count(grepl("^rp", rownames(pbmc.big1)))
# count(!grepl("^mt-", rownames(pbmc.big1)))
gc()
memory.limit(9999999999)

#More plottings; www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html 
as_tibble(
  pbmc.big2[[c("nCount_RNA","nFeature_RNA","percent.mt","Percent.Largest.Gene")]],
  rownames="Cell.Barcode"
) -> qc.metrics2 #change between big1 and big2 (to metrics1 and 2)
#Again, aesthetics in R is not easy to do:
#https://stackoverflow.com/questions/3606697/how-to-set-limits-for-axes-in-ggplot2-r-plots
#https://stackoverflow.com/questions/37950511/r-ggplot2-setting-tick-mark-interval , 
#https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0
# https://felixfan.github.io/ggplot2-remove-grid-background-margin/
qc.metrics1 %>%
  ggplot(aes(Percent.Largest.Gene)) + 
  geom_histogram(binwidth = 0.3, fill="yellow", colour="black") +# ggtitle("Distribution of Percentage Largest Gene") +
  geom_vline(xintercept = 8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  xlab('Percent of Largest Gene (%)')+ylab('Count')+#xlim(0, 30)+
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5), expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))
#8 or 11 with big2 (qc.metrics2), # rm(qc.metrics)

subset(
  pbmc.big1,
  nFeature_RNA>500 & #iteroi kuvaajasta ja yll‰olevalla kaavalla
    nFeature_RNA < 4600 &  #kuten edell‰
    percent.mt <15 & #75% quantile+1sd
    Percent.Largest.Gene < 8 & #vain  yll‰olevalla kaavalla
    percent.ribo >8 &
    #5% quantile-1sd
    percent.hb <0.84  #75% quantile-1sd
) -> pbmc.big1
subset(
  pbmc.big2,
  nFeature_RNA>680 & #iteroi kuvaajasta ja yll‰olevalla kaavalla
    nFeature_RNA < 5000 &  #kuten edell‰
    percent.mt <6 & #75% quantile+1sd
    Percent.Largest.Gene < 11 & #vain  yll‰olevalla kaavalla
    percent.ribo >5 & #5% quantile-1sd
    percent.hb <0.23  #75% quantile
) -> pbmc.big2

#The order is this because you want over mithocondrial cells first away, 
# and then all those genes away from dataset (but not the cells per se) 
dim(pbmc.big1)
#[1] 18297 31738, 18419 19590 (why different dimensions?)
# Filter MALAT1
pbmc.big1 <- pbmc.big1[!grepl("Malat1", rownames(pbmc.big1)), ]
# Filter Mitocondrial
pbmc.big1 <- pbmc.big1[!grepl("^mt-", rownames(pbmc.big1)), ]
# Filter Ribosomal gene (optional if that is a problem on your data) 
pbmc.big1=pbmc.big1[!grepl('^Rps', rownames(pbmc.big1)), ] #remember the equations..
pbmc.big1=pbmc.big1[!grepl('^Rpl', rownames(pbmc.big1)), ]
# Filter Hemoglobin gene
pbmc.big1 <- pbmc.big1[!grepl("^Hb", rownames(pbmc.big1)), ]
dim(pbmc.big1) #head(pbmc.big1@meta.data)
#[1] 18297 31738, 18394 19590 (old..?: remember the subsetting!)

dim(pbmc.big2)
#[1] 19599 27494, 19599 18174 (old?)
# Filter MALAT1
pbmc.big2 <- pbmc.big2[!grepl("Malat1", rownames(pbmc.big2)), ]
# Filter Mitocondrial
pbmc.big2 <- pbmc.big2[!grepl("^mt-", rownames(pbmc.big2)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) 
pbmc.big2=pbmc.big2[!grepl('^Rps', rownames(pbmc.big2)), ] #remember the equations..
pbmc.big2=pbmc.big2[!grepl('^Rpl', rownames(pbmc.big2)), ]
# Filter Hemoglobin gene
pbmc.big2 <- pbmc.big2[!grepl("^Hb", rownames(pbmc.big2)), ]
dim(pbmc.big2)
#[1] 19474 27494, 19573 18174 (old?)

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html
# pbmc1 <- logNormCounts(pbmc.big1)
# dec <- modelGeneVar(pbmc1, block = sce.filt$sample)
# hvgs = getTopHVGs(dec, n = 2000)
# pbmc1 <- runPCA(pbmc1, subset_row = hvgs)
# pbmc1 <- runUMAP(pbmc1, pca = 10)
# # run computeDoubletDensity with 10 principal components.
# pbmc1 <- scDblFinder(pbmc1, dims = 10)
# pbmc.big1 = pbmc.big1[, pbmc1$scDblFinder.score < 2] # This probably not the way to remove doublets..?
#For number of cells per replica or later cluster:
table ( Idents(pbmc.big1) )
table ( Idents(pbmc.big2) )

#Second backup save:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1f_16821.rds")
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1f_16821.rds") #you need later..
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2f_16821.rds")
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2f_16821.rds") #you need later..
gc()
memory.limit(9999999999)
#Earliest stage to calculate the amount of cells with certain gene expressions
#https://github.com/satijalab/seurat/issues/371
#Test Hpl/s or Hbthat you do not have (1.7.21)
# sum(grepl('^Hb', rownames(pbmc.big1)))
# sum(GetAssayData(object = pbmc.big2, slot = "data")["Sox2",]>0)
#Check the subset of cells within your tissue label:
#https://github.com/satijalab/seurat/issues/349
# length(WhichCells(pbmc.big1, slot = 'data', expression = Sox2 > 0 & Pitx2))

#2) Identification of highly variable genes and scaling the data:
#Normalize
# pbmc.big1 <- NormalizeData(pbmc.big1)
#Alternative way to normalize
# NormalizeData(pbmc.big1, normalization.method = "LogNormalize") -> pbmc.big1 #the below is better
NormalizeData(pbmc.big1, normalization.method = "CLR") -> pbmc.big1
NormalizeData(pbmc.big2, normalization.method = "CLR") -> pbmc.big2
#If you can divide between the different types of data
#You may also need cell cycle scoring
# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Cell_Cycle_Scoring
# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Automated_Cell_Type_Annotation
CellCycleScoring(pbmc.big1, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc.big1
as_tibble(pbmc.big1[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
CellCycleScoring(pbmc.big2, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc.big2
as_tibble(pbmc.big2[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

#This 500 was somewhere given as default, I wonder why..2000 should be better
pbmc.big1 <- FindVariableFeatures(pbmc.big1, selection.method = "vst", nfeatures = 2000)
pbmc.big2 <- FindVariableFeatures(pbmc.big2, selection.method = "vst", nfeatures = 2000) 
# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc.big1), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc.big1)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

#We can now access the normalized data in data@assays$RNA@data. 
#We can use this to show that we can get a list of the most highly expressed genes.
apply(pbmc.big1@assays$RNA@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
apply(pbmc.big2@assays$RNA@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
# head(gene.expression, n=50)
#Scaling
#Let's check regressing the cell cycles.. 
#https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html #vars.to.regress = c("S.Score", "G2M.Score")
all.genes <- rownames(pbmc.big1)
pbmc.big1 <- ScaleData(pbmc.big1,  features = all.genes)#vars.to.regress = c("S.Score", "G2M.Score"),
all.genes <- rownames(pbmc.big2)
pbmc.big2 <- ScaleData(pbmc.big2,  features = all.genes)#vars.to.regress = c("S.Score", "G2M.Score"),

#Third backup save:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1b_16821.rds") #check the file name..
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1b_16821.rds") #you need later..
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2b_16821.rds")
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2b_16821.rds") #you need later..
# rm(pbmc.big1,pbmc.big2,pbmc.big1_m)
# gc()
memory.limit(9999999999)
# https://stackoverflow.com/questions/5171593/r-memory-management-cannot-allocate-vector-of-size-n-mb

#3) PCA

#50 is default PC value, if you need more add npcs
#https://www.rdocumentation.org/packages/Seurat/versions/4.0.3/topics/RunPCA
pbmc.big1 <- RunPCA(pbmc.big1, features = VariableFeatures(object = pbmc.big1), npcs=75) #Could be e.g. 75 or 60 if you have less cells..
pbmc.big2 <- RunPCA(pbmc.big2, features = VariableFeatures(object = pbmc.big2), npcs=75)
# Examine and visualize PCA results a few different ways
#print(pbmc.big1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc.big1, dims = 1:2, reduction = "pca")
#DimPlot(pbmc.big1, reduction = "pca")
# DimPlot(pbmc.big1,reduction="pca",label = TRUE)+ggtitle("PC1 vs PC2 with Clusters")
#Fourth backup save:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1pc_16821.rds")
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1pc_16821.rds") #you need later..
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2pc_16821.rds")
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2pc_16821.rds") #you need later..
gc()
memory.limit(9999999999)

#4) Dimensionality of datasets

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc.big1 <- JackStraw(pbmc.big1, num.replicate = 100, dims=70) #Very slow (around 55 min with dim=70), needed? Yes, needed, check below:
pbmc.big1 <- ScoreJackStraw(pbmc.big1, dims = 1:70) #Needed? Less cells, less dimensions
#Plot the data
DimHeatmap(pbmc.big1, dims = 1:15, cells = 500, balanced = TRUE)
JackStrawPlot(pbmc.big1, dims = 1:70) #PC48 has first p-value above 0.05 ->PC47 selected, or 38..
ElbowPlot(object = pbmc.big1, ndims = 70)+ylim(0,8)
# +scale_y_continuous(breaks=seq(0,7,1))
# https://www.biostars.org/p/423306/
pbmc.big2 <- JackStraw(pbmc.big2, num.replicate = 100, dims=70) #Needed? Very slow. Yes, needed, check below:
pbmc.big2 <- ScoreJackStraw(pbmc.big2, dims = 1:70) #Needed?
#Plot the data
DimHeatmap(pbmc.big2, dims = 1:15, cells = 500, balanced = TRUE)
JackStrawPlot(pbmc.big2, dims = 1:70)
ElbowPlot(object = pbmc.big2, ndims = 70)+ylim(0,8)
#Fifth backup save/load:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1JS_16821.rds")
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1JS_16821.rds") #you need later..
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2JS_16821.rds")
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2JS_16821.rds") #you need later..

#5) Clustering of  Cells

#More dimensions are reduced, start with tsne (Non-linear dimensional reduction)
8482 -> saved.seed
set.seed(saved.seed)
  RunTSNE(
    pbmc.big1,
    dims=1:46,
    seed.use = saved.seed, 
    perplexity=200
  ) -> pbmc.big1 #slow, but can do

# DimPlot(pbmc.big1,reduction = "tsne", pt.size = 1) + ggtitle("tSNE with Perplexity 10")
#Non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
#reticulate::py_install(packages ='umap-learn')
#More dimension reductions: UMAP=Uniform Manifold Approximation and Projection
pbmc.big1 <- RunUMAP(pbmc.big1, dims = 1:46) #or 15 etc. check the elbow etc. as above, is this needed if I have already tsne?
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
pbmc.big1 <- FindNeighbors(pbmc.big1, dims = 1:46) #or 10 or 15 or 20, check the elbow etc. as above
#The resolution defines the amount of clusters, higher it is more you have, with 0.5 (default) you should get around 15,
#which sounds biologically more realistic than e.g. 50 with 4
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day3/scRNA_Workshop-PART5.html
new=seq(from = 0.05, to = 5, by = 0.05)
pbmc.big1 <- FindClusters(
  object = pbmc.big1,
  reduction.type = "pca",
  resolution =new, #fast, can do
  dims.use = 1:46,
  save.SNN = TRUE)
head(pbmc.big1$seurat_clusters, n=15)
# DimPlot(pbmc.big1, reduction = "umap")
# DimPlot(pbmc.big1,reduction="tsne",pt.size = 1, label = TRUE, label.size = 7) + ggtitle("tSNE with Perplexity 200")
#Then choose the best resolution
# pbmc.big1 <- FindClusters(pbmc.big1, resolution = 0.5) # unique[]
# head(pbmc.big1, n=3)
# pbmc.big1[, c(13:22)], access to metadata, again not obvious:
#https://satijalab.org/seurat/articles/essential_commands.html
nw= c()
nw2= c()
hei=c()
for (i in 13:112) {#check the dimensions
  nw<-max(as.numeric(unique(pbmc.big1[[]][,i]))) #Thir presumably also for pbmc.big2
  nw2 <-rbind(nw2,nw)}
hei=nw2 #after excel ->0.8 was max, but one below was 0.6, so average is 0.7.. which needs to be calculated:
#Get the right resolution size with group.by as follows:
#https://satijalab.org/seurat/reference/dimplot
pbmc.big1 <- FindClusters(pbmc.big1, resolution = 0.55)
DimPlot(pbmc.big1,reduction="tsne",pt.size = 1, label = TRUE, label.size = 7,group.by = 'RNA_snn_res.0.55') + ggtitle("tSNE with Perplexity 200")
#Sixth backup save and load:
saveRDS(pbmc.big1, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1tsne_16821.rds")

#For the dataset 2 (E14.25), 
8482 -> saved.seed
set.seed(saved.seed)
RunTSNE(
  pbmc.big2,
  dims=1:59,
  seed.use = saved.seed, 
  perplexity=200
) -> pbmc.big2 #slow, but can do
pbmc.big2 <- RunUMAP(pbmc.big2, dims = 1:59) #or 15 etc. check the elbow etc. as above, is this needed if I have already tsne?
# UMAP not needed if tsne, since https://satijalab.org/seurat/articles/essential_commands.html and 
#https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html do not mention both
pbmc.big2 <- FindNeighbors(pbmc.big2, dims = 1:59) #or 10 or 15 or 20, check the elbow etc. as above, new=seq(from = 0.05, to = 5, by = 0.05)
new=seq(from = 0.05, to = 5, by = 0.05) #It is good to have it here in case you need change or just run pbmc.big2 seprately
pbmc.big2 <- FindClusters(
  object = pbmc.big2,
  reduction.type = "pca",
  resolution = new, #fast, can do
  dims.use = 1:59,
  save.SNN = TRUE) #saving below, and loading above
# https://bioinformatics.stackexchange.com/questions/4297/resolution-parameter-in-seurats-findclusters-function-for-larger-cell-numbers
# head(pbmc.big2$seurat_clusters, n=15)
head(pbmc.big2, n=1)
DimPlot(pbmc.big2,reduction="tsne",pt.size = 1, label = TRUE, label.size = 7,group.by = 'RNA_snn_res.0.6') + ggtitle("tSNE with Perplexity 200")#With noncapital letter umap or tsne
nw= c()
nw2= c()
hei2=c()
for (i in 13:112) { #check the dimensions
  nw<-max(as.numeric(unique(pbmc.big2[[]][,i])))
  nw2 <-rbind(nw2,nw)}
hei2=nw2 #after excel ->0.4 is max and 0.2 below so 0.3 is.. or 0.7 as above should be better..
pbmc.big2 <- FindClusters(pbmc.big2, resolution = 0.6) #you need to do this, since the latest is 'seurat_cluster' #This should be 0.7
saveRDS(pbmc.big2, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2tsne_16821.rds")

# Load the data:
# pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1tsne_16821.rds") #you need later..
# pbmc.big1 <- FindClusters(pbmc.big1, resolution = 0.55) # head(pbmc.big1, n=2)
# pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2tsne_16821.rds") #you need later..
# pbmc.big2 <- FindClusters(pbmc.big2, resolution = 0.6) #you need to do this, since the latest is 'seurat_cluster' #This should be 0.7

##This could work also for the PC and the cluster amount determinations and definitions...
# percentage of mitochondrial genes (percent.mito), number of unique molecular identifiers (nUMI), #number of genes expressed (nGene)
# https://chipster.csc.fi/manual/single-cell-seurat-featureplot-v3.html
# https://satijalab.org/seurat/archive/v2.4/pbmc3k_tutorial.html
# https://www.gitmemory.com/issue/NHLBI-BCB/IKAP/2/622987366
# https://github.com/NHLBI-BCB/IKAP
# https://github.com/NHLBI-BCB/IKAP/blob/master/Seurat3_code/example_Seurat3.R
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html
# https://discuss.analyticsvidhya.com/t/how-to-add-a-column-to-a-data-frame-in-r/3278

#6)Plotting the markers and tissues + saving

#Finding differentially expressed features (cluster biomarkers)
# pbmc.big1.markers <- FindAllMarkers(pbmc.big1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(pbmc.big1.markers,file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1.markers_17621.csv")
# https://www.biostars.org/p/410170/
# If the ident.2 parameter is omitted or set to NULL, FindMarkers will test for deferentially expressed features between the group specified by ident.1 and all other cells.
# https://satijalab.org/seurat/archive/v3.0/de_vignette.html
# only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
# Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
# https://www.rdocumentation.org/packages/Seurat/versions/4.0.3/topics/FindMarkers
# #VlnPlot(pbmc.big1, features = c("Sox2", "Pitx2","Foxi3", "Krt14"), slot = "counts", log = TRUE)
# #This does not work before clustering and UMAP procedures ('RunUMAP' mentioned above in this section (7))
# FeaturePlot(pbmc.big1, features = c( "Krt14"))	#"Sox2", "Pitx2","Foxi3",
# #Calculate the amount of genes per cell etc...
# #https://github.com/satijalab/seurat/issues/371
# sum(GetAssayData(object = pbmc.big2, slot = "data")["Sox2",]>0)
# #https://github.com/satijalab/seurat/issues/349
# length(WhichCells(pbmc.big1, slot = 'counts', expression = Pitx2 > 0 ))

#Test this before markers:
table ( Idents(pbmc.big1) )
table ( Idents(pbmc.big2) )
#Markers.. kumpi olikaan oikea luuppi
#Just t‰‰ alempi sill‰ t‰‰ antaa oikeat (ei-kerrannalliset) nimet:
nw= c()
nw2= c()
nww=c()
markers1=c()
# head(pbmc.big1$seurat_clusters, n=15)
# https://stackoverflow.com/questions/33177118/append-a-data-frame-to-a-list
for (i in 0:11) {
  nw<-FindMarkers(pbmc.big1, ident.1 = i, min.pct = 0.07, logfc.threshold = 0.1, test.use = "roc",only.pos = FALSE)
  nw$cluster=i
  nww=list(data.frame(nw))
  nw2 <-list.append(nw2,nww)}
markers1=nw2
# nw= c()
# nw2= c()
# markers1=c()
# # head(pbmc.big1$seurat_clusters, n=15)
# for (i in 0:11) {
#   nw<-FindMarkers(pbmc.big1, ident.1 = i, min.pct = 0.07, logfc.threshold = 0.1, test.use = "roc",only.pos = FALSE)
#   nw$cluster=i
#   nw2 <-rbind(nw2,nw)}
# markers1=nw2 #markers1=markers11, rm(markers2), row.names(data.frame(markers1[7]))
nw2= c()
nww=c()
for (i in 1:12) {
  nww=list(row.names(data.frame(markers1[i])))
  nw2 <-list.append(nw2,nww)}
m1_names=unlist(list(nw2))
  # list(row.names(data.frame(markers1))) #nw2
m1=data.frame()
for (i in 1:12) {
  m1=rbind(m1,data.frame(markers1[i]))}
markers1=m1
# write.csv(markers1,file = "D:/Results/Markers/markers_E11.5_16821.csv") #getwd() # head(markers1[3])
# https://stackoverflow.com/questions/33177118/append-a-data-frame-to-a-list
# write.csv(markers1,file = "D:/Results/Markers/markers_E11.5na.csv") #getwd()
# nw= c()
# nw2= c()
# # head(pbmc.big2$seurat_clusters, n=15)
# markers2=c()
# for (i in 0:17) {
#   nw<-FindMarkers(pbmc.big2, ident.1 = i, min.pct = 0.07, logfc.threshold = 0.1, test.use = "roc",only.pos = FALSE)
#   nw$cluster=i
#   nw2 <-rbind(nw2,nw)}
# markers2=nw2

nw= c()
nw2= c()
nww=c()
markers2=c()
# head(pbmc.big2$seurat_clusters, n=15)
for (i in 0:17) {
  nw<-FindMarkers(pbmc.big2, ident.1 = i, min.pct = 0.07, logfc.threshold = 0.1, test.use = "roc",only.pos = FALSE)
  nw$cluster=i
  nww=list(data.frame(nw))
  nw2 <-list.append(nw2,nww)}
markers2=nw2

nw2= c()
nww=c()
for (i in 1:18) {
  nww=list(row.names(data.frame(markers2[i])))
  nw2 <-list.append(nw2,nww)}
# m2_names=list(row.names(data.frame(markers2))) #nw2
m2_names=unlist(list(nw2))

m2=data.frame()
for (i in 1:18) {
  m2=rbind(m2,data.frame(markers2[i]))}
markers2=m2
# write.csv(markers2,file = "D:/Results/Markers/markers_E14.25_16821.csv")
# m1$X=unlist(m1_names)
# samp1 <- m1[,1:7]
# # rownames(samp1) <- m1[,8]
# markers1=samp1

############# You need add gene names to the previous, sound simple but it is needed for sub settings:
# https://stackoverflow.com/questions/5555408/convert-the-values-in-a-column-into-row-names-in-an-existing-data-frame/55390036
# markers1[1:5,1:7]
#Meniskˆh‰n t‰‰ n‰in..
# https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists
markers1 <- cbind(markers1, markers1$pct.1*2*(markers1$myAUC-0.5)+(1-markers1$pct.2)/3)
colnames(markers1)[8] <- "Positive Markers"
markers1 <- cbind(markers1, (1-markers1$pct.1)*2*(0.5-markers1$myAUC)+markers1$pct.2/4)
colnames(markers1)[9] <- "Negative Markers"
markers2 <- cbind(markers2, markers2$pct.1*2*(markers2$myAUC-0.5)+(1-markers2$pct.2)/3)
colnames(markers2)[8] <- "Positive Markers"
markers2 <- cbind(markers2, (1-markers2$pct.1)*2*(0.5-markers2$myAUC)+markers2$pct.2/4)
colnames(markers2)[9] <- "Negative Markers"

markers1 <- cbind(markers1, markers1$pct.1-markers1$pct.2)
colnames(markers1)[10] <- "Cluster Difference"
markers2 <- cbind(markers2, markers2$pct.1-markers2$pct.2)
colnames(markers2)[10] <- "Cluster Difference"

markers1 <- cbind(markers1, (markers1$pct.1-markers1$pct.2)*markers1$myAUC)
colnames(markers1)[11] <- "Goodness of a Cluster Marker"
markers2 <- cbind(markers2, (markers2$pct.1-markers2$pct.2)*markers2$myAUC)
colnames(markers2)[11] <- "Goodness of a Cluster Marker"
# markers1 <- markers1[,c(1,2,3,4,5,6,8,9,7)] #Switching columns
# markers1 <- markers1[,-c(13,14,15)] #In case of too many columns you can remove one like this


#The first (E11.5)
# https://www.statmethods.net/management/operators.html
# markers1 <- markers1[,-c(12,13)]
# https://statisticsglobe.com/r-warning-message-condition-length-only-first-element-used
vec=c()
for(i in 1:length(markers1$myAUC)) { 
  if (markers1$myAUC[i] >=0.5) {
    if (markers1$pct.1[i]-markers1$pct.2[i] >= 0.2) {
      vec[i] <- markers1$pct.1[i]*(markers1$myAUC[i]-0.5)+0.5+(markers1$pct.1[i]-markers1$pct.2[i])/3

    } else if (markers1$pct.1[i]-markers1$pct.2[i] < 0.2) {
    vec[i]  <- markers1$pct.1[i]*(markers1$myAUC[i]-0.5)+0.5+abs(markers1$pct.1[i]-markers1$pct.2[i])/5}

} else if (markers1$myAUC[i] < 0.5) {
  vec[i] <-  0
} else {
  vec[i] <-  0}}
markers1 <- cbind(markers1,vec)
colnames(markers1)[12] <- "Ultimate Positive"
vec=c()
for(i in 1:length(markers1$myAUC)) { 
  if (markers1$myAUC[i] < 0.5) {
    if (markers1$pct.2[i]-markers1$pct.1[i] >= 0.2) {
      vec[i] <- (1-markers1$pct.1[i])*(0.5-markers1$myAUC[i])+0.5+(markers1$pct.2[i]-markers1$pct.1[i])/3
      
    } else if (markers1$pct.2[i]-markers1$pct.1[i] < 0.2) {
      vec[i]  <- (1-markers1$pct.1[i])*(0.5-markers1$myAUC[i])+0.5+abs(markers1$pct.2[i]-markers1$pct.1[i])/5}
    
  } else if (markers1$myAUC[i]  >= 0.5) {
    vec[i] <-  0
  } else {
    vec[i] <-  0}}
markers1 <- cbind(markers1,vec)
colnames(markers1)[13] <- "Ultimate Negative"

#For the other (E14.25)
for(i in 1:length(markers2$myAUC)) { 
  if (markers2$myAUC[i] >=0.5) {
    if (markers2$pct.1[i]-markers2$pct.2[i] >= 0.2) {
      vec[i] <- markers2$pct.1[i]*(markers2$myAUC[i]-0.5)+0.5+(markers2$pct.1[i]-markers2$pct.2[i])/3
      
    } else if (markers2$pct.1[i]-markers2$pct.2[i] < 0.2) {
      vec[i]  <- markers2$pct.1[i]*(markers2$myAUC[i]-0.5)+0.5+abs(markers2$pct.1[i]-markers2$pct.2[i])/5}
    
  } else if (markers2$myAUC[i] < 0.5) {
    vec[i] <-  0
  } else {
    vec[i] <-  0}}
markers2 <- cbind(markers2,vec)
colnames(markers2)[12] <- "Ultimate Positive"
vec=c()
for(i in 1:length(markers2$myAUC)) { 
  if (markers2$myAUC[i] < 0.5) {
    if (markers2$pct.2[i]-markers2$pct.1[i] >= 0.2) {
      vec[i] <- (1-markers2$pct.1[i])*(0.5-markers2$myAUC[i])+0.5+(markers2$pct.2[i]-markers2$pct.1[i])/3
      
    } else if (markers2$pct.2[i]-markers2$pct.1[i] < 0.2) {
      vec[i]  <- (1-markers2$pct.1[i])*(0.5-markers2$myAUC[i])+0.5+abs(markers2$pct.2[i]-markers2$pct.1[i])/5}
    
  } else if (markers2$myAUC[i]  >= 0.5) {
    vec[i] <-  0
  } else {
    vec[i] <-  0}}
markers2 <- cbind(markers2,vec)
colnames(markers2)[13] <- "Ultimate Negative"

# Riit‰iskˆh‰n n‰‰?..
# write.csv(row.names(data.frame(markers1[7])),file = "D:/Results/Markers/m1_6g.csv")
# data.frame(m1_names[12])
markers1$X=unlist(m1_names) #this works..
markers2$X=unlist(m2_names) #this works.. :)
markers1=markers1[,c(14,1,2,3,4,5,6,7,8,9,10,11,12,13)]
markers2=markers2[,c(14,1,2,3,4,5,6,7,8,9,10,11,12,13)]
# markers1[,'X.1']=markers1[,'X']
names(markers1)[1]='Genes'
# markers1=markers1[,1:14]
# markers2[,'X.1']=markers2[,'X']
names(markers2)[1]='Genes'
# markers2=markers2[,1:14]
# m11=markers1
# m22=markers2
write.csv(markers1,file = "D:/Results/Markers/markers_E11.5_16821.csv")
write.csv(markers2,file = "D:/Results/Markers/markers_E14.25_16821.csv")
#For panther, c7 in one (E11.5), or cluster2 in two (E14.25)
# write.csv(markers1,file = "D:/Results/Markers/markers_1_c7.csv")
# m110=markers1["avg_log2FC"] #^2 #note the correct dataframe..
# write.table(m110,file = "D:/Results/Genes/c110_E115.txt",row.names = TRUE, col.names = FALSE,quote = FALSE, sep = '\t') 
# For cluster 2 in two (E14.25); Panther Analysis details:
# Mapped IDs:	1686 out of 1718
# Unmapped IDs:	68
# Multiple mapping information:	54, ok. unless otherwise stated.. rm()

#For later purposes:
# rm(markers1,markers2)
# markers1=read.csv("D:/Results/Markers/markers_E11.5_16821.csv", header = TRUE)
# markers2=read.csv("D:/Results/Markers/markers_E14.25_16821.csv", header = TRUE)


# samp1 <- markers1[,-1]
# rownames(samp1) <- markers1[,'X']
# markers1=samp1
# samp2 <- markers2[,-1]
# rownames(samp2) <- markers2[,1]
# markers2=samp2
# rm(samp1,samp2)
# dim(data.frame(m1_names[12]))[1]
# length(row.names(data.frame(filter(markers1, cluster == 11))))

#This works for (E11.5)
# dir.create("D:/Results/Markers/Tests")
# https://dplyr.tidyverse.org/reference/filter.html
# filter(markers1, cluster == 5)[order(filter(markers1, cluster == 5)[,12], decreasing=TRUE),]
#This works, list.append and list()
# https://www.programmingr.com/fast-r-append-list/
# https://www.geeksforgeeks.org/list-of-dataframes-in-r/
myList<-list()
for(i in 0:11) {
  myList<-list.append(myList,data.frame(filter(markers1, cluster == i)[order(filter(markers1, cluster == i)[,13], decreasing=TRUE),]))}
# https://statisticsglobe.com/split-vector-into-chunks-in-r
# myList2=split(myList,             # Applying split() function
#               cut(seq_along(myList),
#                   12,
#                   labels = TRUE)) #not needed if you do list.append!
# https://statisticsglobe.com/r-write-read-multiple-csv-files-for-loop
# install.packages("xlsx")
library("xlsx")
for(i in 1:12) {                              # Head of for-loop
  write.xlsx(data.frame(myList[i]),                              # Write CSV files to folder
             paste("D:/Results/Markers/Tests/",i-1,".xlsx"),
             row.names = FALSE)} #can be also with row names: row.names = TRUE

#This works for (E14.25)
# dir.create("D:/Results/Markers/Tests2")
# https://dplyr.tidyverse.org/reference/filter.html
myList<-list()
for(i in 0:17) {
  myList<-list.append(myList,data.frame(filter(markers2, cluster == i)[order(filter(markers2, cluster == i)[,13], decreasing=TRUE),]))}
for(i in 1:18) {                              # Head of for-loop
  write.xlsx(data.frame(myList[i]),                              # Write CSV files to folder
             paste("D:/Results/Markers/Tests2/",i-1,".xlsx"),
             row.names = FALSE)} #can be also with row names: row.names = TRUE
#Then do this:
# https://www.ablebits.com/office-addins-blog/2017/11/08/merge-multiple-excel-files-into-one/
# https://www.rdocumentation.org/packages/Seurat/versions/4.0.3/topics/FindMarkers
# https://discuss.analyticsvidhya.com/t/how-to-add-a-column-to-a-data-frame-in-r/3278

#How to check if cluster has certain gene
# FeaturePlot(pbmc.big1, features = c("Pitx2"))
# VlnPlot(pbmc.big1,features=c("Sox2"))
# VlnPlot(pbmc.big1,features=c("Krt14"))
#For number of cells per cluster:
table ( Idents(pbmc.big1) )
table ( Idents(pbmc.big2) )
#Cell Cycles
table (Idents(pbmc.big1),pbmc.big1$Phase)
table (Idents(pbmc.big2),pbmc.big2$Phase)
# VlnPlot(pbmc.big1, features = c("Phase"))
# https://github.com/satijalab/seurat/issues/354

# krts=rownames(avgEt2[grep('^Krt', rownames(avgEt2)),])
#How to get all the groups that have at least Krts and Epcam. One can test to add e.g. Pitx2 and Sox2 and Foxi3
#This is for pbmc.big1  
# krts=rownames(avgEt2[grep('^Krt', rownames(avgEt2)),])
avgEt=AverageExpression(object = pbmc.big1)
# https://stackoverflow.com/questions/43425300/list-of-lists-to-matrix
avgEt2 <- do.call("cbind",avgEt) #This is needed!
krts2=rownames(avgEt2[grep('^Krt', rownames(avgEt2)),])
krts22=rownames(avgEt2[grep('^Wnt', rownames(avgEt2)),])
# krts2=append(krts2, c("Epcam","Pitx2","Sox2","Foxi3")) #Sostdc1.. and more:
krts2=append(krts2, c("Epcam","Pitx2","Sox2","Foxi3","Sostdc1", "Dlx2", "Lef1", "p63", "Shh", "Noggin", "Sema3f", 
                      "Fgf8", "Hh", "Axin2", "FzD6", "Sp6", "Edar","TGFB1", "ACTA2","Lgr5", "Dkk3", "Igfp5", 
                      "Bmi1", "ABCG2", "Oct3/4", "Yap", "Gli1", "Lrig1", "Fgf10", "p21", "Dkk4", "Fgf9", "Fgf20",
                      'Ptma', 'Tpt1', 'Enam', 'Amelx','Mmp20', 'Amtn', 'Klk4', 'Dbi', 'Acta','Actb','Tmsb4x','Fth1','Pttg1', 'Atf3', 'Cldn10'))
krts2=append(krts2,krts22)
# https://hbctraining.github.io/Intro-to-R/lessons/06_matching_reordering.html
krts2=krts2[krts2 %in% rownames(avgEt2)]
# sum(GetAssayData(object = pbmc.big1, slot = "data")["Krt10",]>0) #or14
# https://stackoverflow.com/questions/18933187/how-to-select-some-rows-with-specific-rownames-from-a-dataframe
krt_DE1=markers1[krts2,] #Check markers1 to 2 and vice versa
kde1=sort(unique(krt_DE1$cluster))
# https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
krt_DE1 = na.omit(krt_DE1)
krt_DE1

avgEt=AverageExpression(object = pbmc.big2)
avgEt2 <- do.call("cbind",avgEt) #This is needed!
krts2=rownames(avgEt2[grep('^Krt', rownames(avgEt2)),])
krts22=rownames(avgEt2[grep('^Wnt', rownames(avgEt2)),])
krts2=append(krts2, c("Epcam","Pitx2","Sox2","Foxi3","Sostdc1", "Dlx2", "Lef1", "p63", "Shh", "Noggin", "Sema3f",
                      "Fgf8", "Hh", "Axin2", "FzD6", "Sp6", "Edar","TGFB1", "ACTA2","Lgr5", "Dkk3", "Igfp5",
                      "Bmi1", "ABCG2", "Oct3/4", "Yap", "Gli1", "Lrig1", "Fgf10", "p21", "Dkk4", "Fgf9", "Fgf20",
                      'Ptma', 'Tpt1', 'Enam', 'Amelx','Mmp20', 'Amtn', 'Klk4', 'Dbi', 'Acta','Actb','Tmsb4x','Fth1','Pttg1', 'Atf3', 'Cldn10'))
#PTMA and TPT1, should be here as well, and Dbi
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=TPT1
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=PTMA
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=DBI&keywords=Dbi
#and ameoblasts:
# ameloblast (Enam, Amelx, Mmp20, Amtn, Klk4) 
# Mesenchyme, Col1A
# https://pubmed.ncbi.nlm.nih.gov/29405385/
#Blood epithelia, (Actb), and Tmsb4x and Fth13 
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=TMSB4X&keywords=Tmsb4x
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=ACTB&keywords=Actb
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=FTH1&keywords=Fth1
krts2=append(krts2,krts22)
krts2=krts2[krts2 %in% rownames(avgEt2)]
krt_DE2=markers2[krts2,] #Check markers1 to 2 and vice versa
kde2=sort(unique(krt_DE2$cluster))#Check this and kde1
krt_DE2 = na.omit(krt_DE2)
krt_DE2
#If you need print more than around 20 markers
# https://stackoverflow.com/questions/23188900/display-print-all-rows-of-a-tibble-tbl-df
# krt_DE2 %>% tibble::as_tibble() %>% print(n=Inf)
# avgEt=AverageExpression(object = pbmc.big1)
# https://stackoverflow.com/questions/43425300/list-of-lists-to-matrix
# avgEt2 <- do.call("cbind",avgEt) #This is needed!
# https://dzone.com/articles/learn-r-how-extract-rows
# Epit=avgEt2[krts2,] #for second

#Let this be the lower limit for epithelia:
# Epit_Exp_Min=quantile(Epit[ !Epit==0 ], p=0.1)
# Epit_Exp_Min=quantile(Epit[ !Epit==0 ], p=0.5) #this is good or check histo
# hist( Epit[ !Epit==0 ], breaks=2000,xlim=c(0,0.1))
# https://www.datacamp.com/community/tutorials/make-histogram-basic-r
# krt_clusters=avgEt2[krts2,]>Epit_Exp_Min
# length(krt_clusters[,1])
# sum(krt_clusters[,1], na.rm = TRUE) 

# https://www.r-bloggers.com/2015/12/how-to-write-the-first-for-loop-in-r/
for (i in 1:12){print(sum(krt_clusters[,i], na.rm = TRUE))} #check the cluster size
# nk.raw.data <- as.matrix(GetAssayData(nrd)[, WhichCells(nrd, ident = "1")])
# length(nk.raw.data)
# sum(GetAssayData(object = pbmc.big1, slot = "data", Idents(1))["Sox2",]>0)
#for both
table ( Idents(pbmc.big1) ) 
table ( Idents(pbmc.big2) ) 

#This works the best..
#https://www.r-bloggers.com/2012/10/error-handling-in-r/
nrd=subset(pbmc.big2, subset = Vim > 0) #If you need more... develop a strategy
b= c()
nw2=c()
for(i in 0:17) {
  b=tryCatch( print( dim(subset(x = nrd, idents = i))[2] ) ,
            warning = function(w) {print(paste("No cells in cluster", i)); 
              log(-i)},
            error = function(e) {print(0);
              0}) 
  nw2<-append(nw2,list(b))}
nw2=unlist(nw2)
nw2

### Testing an alternative way to analyse markers, MGFR: marker gene finder for rna-seq
# https://bioconductor.org/packages/release/bioc/html/MGFR.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("MGFR")
# library("MGFR")
# ls("package:MGFR")
# data("ref.mat")
# dim(ref.mat)
# colnames(ref.mat)
# markers.list <- getMarkerGenes.rnaseq(ref.mat, class.vec = colnames(ref.mat),
#                                      samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
# names(markers.list)
# # show the first 20 markers of liver
# markers.list[["liver_markers"]][1:20]
# sessionInfo()



#Monocle test for pseudotime analysis
# https://github.com/satijalab/seurat-wrappers 
# https://satijalab.org/signac/articles/monocle.html
# http://cole-trapnell-lab.github.io/monocle-release/docs/#trajectory-step-2-reduce-data-dimensionality
#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# rm(pbmc.big1,pbmc.big2,pbmc.big1_m)
pbmc.big1 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big1tsne_16821.rds") #you need later..
# pbmc.big1 <- FindClusters(pbmc.big1, resolution = 0.55)
pbmc.big2 <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/pbmc.big2tsne_16821.rds") #you need later..
# pbmc.big2 <- FindClusters(pbmc.big2, resolution = 0.6) #you need to do this, since the latest is 'seurat_cluster' #This should be 0.7

#This is ok 9721 Tikka:
# https://github.com/cole-trapnell-lab/monocle3/issues/262
seurat_object=c()
seurat_object <- identity(pbmc.big1)
gene_symbol <- as.list(org.Mm.egSYMBOL)
raw_count_data <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")
class(raw_count_data)
cells_info <- seurat_object@meta.data
gene_name <- gene_symbol[rownames(raw_count_data)]  
gene_name <- sapply(gene_name, function(x) x[[1]][1])
#preparing cds
gene_name <- ifelse(is.na(gene_name), names(gene_name), gene_name)
gene_short_name <- gene_name
gene_id <- rownames(raw_count_data)
genes_info <- cbind(gene_id, gene_short_name)
genes_info <- as.data.frame(genes_info)
rownames(genes_info) <- rownames(raw_count_data)
cds=c()
cds <- new_cell_data_set(expression_data = raw_count_data,
                         cell_metadata = cells_info,
                         gene_metadata = genes_info)
# cds <- as.cell_data_set(pbmc.big1) #n‰in tehtyn‰ pseudoaika-plotti n‰ytt‰‰ erilaiselta kuin aiempi cds..

cds <- preprocess_cds(cds, num_dim = 46)
cds <- reduce_dimension(cds)
# cds <- reduce_dimension(cds,reduction_method = 'tSNE')#head(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = TRUE)
cds <- order_cells(cds)#Select Mesenchymal cells points in the middle and epithelial otherwise (in the surface)
# setwd(dir='D:/Data for Seurat Analysis/')
save(cds,file="cds.Rdata")
# load(file="cds.Rdata")
# plot_cells(cds, reduction_method="tSNE",label_cell_groups=TRUE)+ theme(legend.position = "right")
plot_cells(cds,          #or cds or pbmc.big1_m
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) #Check this before gene_pseudotime analysis

# https://davetang.org/muse/2017/10/01/getting-started-monocle/ , e.g. fData(cds)$gene_id
ciliated_cds_pr_test_res=c()
# trace('calculateLW', edit = T, where = asNamespace("monocle3")) line93: tmp <- rbind(tmp, cur_tmp)
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
#save this variable, it takes 27min to perform this 
save(ciliated_cds_pr_test_res,file="ciliated_cds_pr_test_res.Rdata")
# load(file="ciliated_cds_pr_test_res.Rdata")

# https://groups.google.com/g/monocle-3-users/c/tBjYuAxwyEo
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
# AFD_genes <- c("Pitx2","Cdh5","Ccl12","Smpd3","Chchd10","Dsp","Fcer1g",'Fn1') #"Krt14","Epcam",
#IEE/OEE, Endothelial Cells, Leukocytes, Odontoblasts, SI/SR, Ameloblast (Dsp), Blood, Mesenchyme ("Fn1" or Vim)
AFD_genes1 <- c("Sox2") #"Epcam","Pitx2","Krt14", "Igfbp5","Sox2"
AFD_lineage_cds <- cds[rowData(cds)$gene_id %in% AFD_genes1,colData(cds)$seurat_clusters %in% c(10)] #c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,11)
#check genes!
# colData(cds)$seurat_clusters %in% 10] #gene_short_name away..
# In "plot_genes_in_pseudotime" min_expression is crucial for seeing anything, 
#it needs to be low for many important genes such as Sox2: 0.0001, or even 0.000003.., or 1e-16
# Could be also color_cells_by="embryo.time.bin", or G/M/S score..# +ylim(0,0.015) #+ggtitle("Sox2")
# head(pData(AFD_lineage_cds)) #the column for G/S score is 'phase': color_cells_by="Phase" works
plot_genes_in_pseudotime(AFD_lineage_cds, color_cells_by="pseudotime", min_expr=0.001,label_by_short_name=FALSE) #Here is the plot
+ylim(0,1)
# +geom_point()+
# https://cole-trapnell-lab.github.io/monocle3/docs/differential/
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c( 0.01 ))#0,10^seq(-6,-1), check later..
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$RNA_snn_res.0.55)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
out=pheatmap::pheatmap(agg_mat,
                       scale="column", clustering_method="ward.D2")
# https://www.biostars.org/p/287512/
# out$tree_row[["order"]]#jehei #length(out$tree_row[["order"]])
# out$tree_row[["order"]][66-9]
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(19,17)), #check..  7, 30.. 19, 29
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
mat=gene_module_df %>% filter(module %in% c(19,29))
mat[order(mat[,5], decreasing = TRUE),]
# 1 C1ql3  19     1            3.91  5.11
# 2 Pdgfa  19     1            3.93  5.08
mat[order(mat[,4], decreasing = TRUE),]
# 1 Lbx1    19     1            4.06  4.47
# 2 Hoxb3   19     1            4.05  4.50
#replace monocle dimensions with seurat.. needed?
# seurat_object <- RunTSNE(seurat_object, dims = 1:46)
cds@reducedDims$tSNE <-  seurat_object@reductions$tsne@cell.embeddings
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122333/

#preparing cds
seurat_object2=c()
seurat_object2 <- identity(pbmc.big2)
gene_symbol <- as.list(org.Mm.egSYMBOL)
raw_count_data <- GetAssayData(seurat_object2, assay = "RNA", slot = "counts")
class(raw_count_data)
cells_info <- seurat_object2@meta.data
gene_name <- gene_symbol[rownames(raw_count_data)]  
gene_name <- sapply(gene_name, function(x) x[[1]][1])
gene_name <- ifelse(is.na(gene_name), names(gene_name), gene_name)
gene_short_name <- gene_name
gene_id <- rownames(raw_count_data)
genes_info <- cbind(gene_id, gene_short_name)
genes_info <- as.data.frame(genes_info)
rownames(genes_info) <- rownames(raw_count_data)
cds2=c()
cds2 <- new_cell_data_set(expression_data = raw_count_data,
                         cell_metadata = cells_info,
                         gene_metadata = genes_info)
# cds2 <- as.cell_data_set(pbmc.big2)
cds2 <- preprocess_cds(cds2, num_dim = 59)
cds2 <- reduce_dimension(cds2) # cds2 <- reduce_dimension(cds2,reduction_method = 'tSNE')# head(cds2), this tSNE is not working here..
cds2 <- cluster_cells(cds2)
cds2 <- learn_graph(cds2, use_partition = TRUE)
cds2 <- order_cells(cds2)
save(cds2,file="cds2.Rdata")
# load(file="cds2.Rdata")
plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

#Change Matrix::rBind to just rbind with/in (entering): trace('calculateLW', edit = T, where = asNamespace("monocle3"))
ciliated_cds_pr_test_res2 <- graph_test(cds2, neighbor_graph="principal_graph", cores=4) #If you have a pseudotime) image, do not (!) enlarge it during calculation,
#It crashes the session.. and check cores.. 4  to 1, save the variable, it takes 27min
save(ciliated_cds_pr_test_res2,file="ciliated_cds_pr_test_res2.Rdata")
# load(file="ciliated_cds_pr_test_res2.Rdata")
pr_deg_ids2 <- row.names(subset(ciliated_cds_pr_test_res2, q_value < 0.05))
# AFD_genes <- c("Krt14","Epcam","Cdh5","Ccl12","Smpd3","Chchd10","Dsp","Fcer1g","Vim")
#IEE/OEE, Endothelial Cells, Leukocytes, Odontoblasts, SI/SR, Ameloblast (Dsp), Blood, Mesenchyme
AFD_genes2 <- c("Sox2") #"Epcam","Pitx2","Krt14", "Igfbp5","Sox2"
AFD_lineage_cds2 <- cds2[row.names(ciliated_cds_pr_test_res2) %in% AFD_genes2,
                         colData(cds2)$seurat_clusters %in% c(2)] 
#check genes! c(0, 1, 3, 4, 5, 6, 7, 8, 9,10,11,12, 13, 14, 15, 16, 17)
# ] #gene_short_name away..
# In "plot_genes_in_pseudotime" min_expre is crusial, it needs to be low for many important genes such as Sox2: 0.02, check this
# Could be also color_cells_by="embryo.time.bin", or G/M/S score..# +ylim(0,0.015) #+ggtitle("Sox2")
# head(pData(AFD_lineage_cds)) #the column for G/S score is 'phase': color_cells_by="Phase" works
plot_genes_in_pseudotime(AFD_lineage_cds2, color_cells_by="pseudotime", min_expr=0.015,label_by_short_name=FALSE)+ylim(0.1,0.4)

gene_module_df2 <- find_gene_modules(cds2[pr_deg_ids2,], resolution=c( 0.01 ))#0,10^seq(-6,-1), check later..
cell_group_df2 <- tibble::tibble(cell=row.names(colData(cds2)), 
                                cell_group=colData(cds2)$RNA_snn_res.0.6)
agg_mat2 <- aggregate_gene_expression(cds2, gene_module_df2, cell_group_df2)
row.names(agg_mat2) <- stringr::str_c("Module ", row.names(agg_mat2))
Out=pheatmap::pheatmap(agg_mat2,
                   scale="column", clustering_method="ward.D2")
# https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-August-Advanced-scRNAseq/master/data_analysis/adv_scrnaseq_monocle.Rmd
plot_cells(cds2,
           genes=gene_module_df2 %>% filter(module %in% c(1,5)), #check..  7, 30..
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
# gene_module_df2 %>% filter(module %in% c(27,38))
# (gene_module_df2 %>% filter(module %in% c(27,38)))$id
mat2=gene_module_df2 %>% filter(module %in% c(1,5))
mat2[order(mat2[,4], decreasing = TRUE),] #Was there minus values?
mat2[order(mat2[,5], decreasing = TRUE),]

#More about markers.. comparison
marker_test_res1 <- top_markers(cds, group_cells_by="seurat_clusters", #can this be cds?
                                reference_cells=1000, cores=8) #could be also pbmc.big1_m or other
marker_test_res2 <- top_markers(cds2, group_cells_by="seurat_clusters", #can this be cds?
                                reference_cells=1000, cores=8) #could be also pbmc.big1_m or other

write.csv(marker_test_res1,file = "D:/Results/Markers/markers_E11.5_17821_alt.csv")
write.csv(marker_test_res2,file = "D:/Results/Markers/markers_E14.25_17821_alt.csv")

# marker_test_res1=read.csv("D:/Results/Markers/markers_E11.5_17821_alt.csv", header = TRUE)
# marker_test_res2=read.csv("D:/Results/Markers/markers_E14.25_17821_alt.csv", header = TRUE)


# colData(cds)
top_specific_markers <- marker_test_res2 %>%
  filter(fraction_expressing >= 0.1) %>%
  group_by(cell_group) %>%
  top_n(15, pseudo_R2) #If you want to see Sox2 this should be 15 in cds2
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id)) 
'Sox2' %in% top_specific_marker_ids #which('Sox2' == top_specific_marker_ids)
# https://stackoverflow.com/questions/1169248/test-if-a-vector-contains-a-given-element
rowData(cds)$gene_short_name <- row.names(rowData(cds))
plot_genes_by_group(cds,top_specific_marker_ids,group_cells_by="seurat_clusters",
                    ordering_type="maximal_on_diag",max.size=1.6)+ theme(text = element_text(size = 6.5))  
# https://cole-trapnell-lab.github.io/monocle3/docs/clustering/

# Require that markers have at least JS specificty score > 0.5 and
# be significant in the logistic test for identifying their cell type:
garnett_markers <- marker_test_res2 %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)
# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% 
  group_by(gene_short_name) %>%
  filter(n() == 1) #from: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
generate_garnett_marker_file(garnett_markers, file="./marker_file2.txt")

# cds_subset <- preprocess_cds(pbmc.big1_m, num_dim = 50)
# cds_subset = reduce_dimension(cds_subset, max_components = 2)
# cds_subset = cluster_cells(cds_subset)
# # modulated_genes <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 4)
# cds_subset <- learn_graph(cds_subset, use_partition = TRUE)
# plot_cells(cds_subset,color_cells_by = "partition") #
# cds_subset <- order_cells(cds_subset)
plot_genes_in_pseudotime(cds[1:100,])#..

#Other good info for graph_test:
# https://groups.google.com/g/monocle-3-users/c/tBjYuAxwyEo
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# https://github.com/cole-trapnell-lab/monocle-release/issues/295
# https://www.biostars.org/p/493489/

#Enrichment
#GO enrichment:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7490294/#footnote3 #this is the original interest
# http://pantherdb.org/index.jsp #this particularly
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323939/ #info regarding panther
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037419/
genes=read.csv("D:/Codes and Instructions/ll.csv", header = FALSE) #where did you get ll.csv?
# https://intellipaat.com/community/15241/how-to-remove-last-n-characters-from-every-element-in-the-r-vector
# https://stringr.tidyverse.org/reference/str_replace.html #could be ok for later purposes as matching..
# https://www.programmingr.com/examples/neat-tricks/r-match-and-in-operator/ #similarly
# https://statisticsglobe.com/r-extract-first-or-last-n-characters-from-string #sim.
# https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html #sim
genes$V1 = substr(genes$V1,1,nchar(genes$V1)-1)
genes_o=read.csv("D:/Codes and Instructions/l2.txt", header = FALSE)
# str_remove(genes_o, genes)
# match(genes, genes_o)
#https://stackoverflow.com/questions/11134812/how-to-find-the-length-of-a-string-in-r #, which enables:
nchar(substr(genes_o[1,1], 1, 3))   
#Eureka!
# https://stackoverflow.com/questions/3492379/data-frame-rows-to-a-list
x <- as.list(as.data.frame(t(genes_o)))# in the example: c("thing", "something", "some", "else")
test <- as.list(as.data.frame(t(genes)))# in the example: c("thing", "some")
pattern <- paste(test, collapse = "|") # create regex pattern
df1 <- data.frame(genes_o,stringsAsFactors=FALSE)
df2 <- data.frame(genes,stringsAsFactors=FALSE)
#This is close enough:
df1[grepl(pattern, x),] <- substr(df1[grepl(pattern, x),],1,nchar(df1[grepl(pattern, x),])-1)
write.csv(df1,file = "D:/Results/Markers/c2_genes.csv")
write.xlsx(df1,                              # Write CSV files to folder
           paste("D:/Results/Genes/c2_E14.25.xlsx"),
           row.names = TRUE)

m1=markers1["avg_log2FC"]^2
# write.csv(m1,file = "D:/Results/Markers/all1_genesi.txt",row.names = TRUE)
#This eventually worked, but take apostrophes ("s) away (replace with 'nothing')
write.table(m1,file = "D:/Results/Markers/all1_genesa.txt",row.names = TRUE, col.names = FALSE,quote = FALSE, sep = '\t')
#One does need to do this cluster by cluster in loop with all data
#Let's start with the important cluster 2 in data2:
head(markers2[markers2$cluster==2,], n=10)
m2=markers2[markers2$cluster==2,]
m22=m2["avg_log2FC"]^2 #note the correct dataframe..
write.table(m22,file = "D:/Results/Genes/c2_E1425.txt",row.names = TRUE, col.names = FALSE,quote = FALSE, sep = '\t')
genes=read.csv("D:/Results/Genes/c2_c.txt", header = FALSE)
x <- as.list(as.data.frame(t(rownames(m2))))# in the example: c("thing", "something", "some", "else")
test <- as.list(as.data.frame(t(genes)))# in the example: c("thing", "some")
pattern <- paste(test, collapse = "|") # create regex pattern
df1 <- data.frame(x,stringsAsFactors=FALSE)
df2 <- data.frame(test,stringsAsFactors=FALSE)
#This is close enough:
df1[grepl(pattern, x),] <- substr(df1[grepl(pattern, x),],1,nchar(df1[grepl(pattern, x),])-1)
write.xlsx(df1,                              # Write CSV files to folder
           paste("D:/Results/Genes/c2_E1425ok2.xlsx"),
           row.names = TRUE)
# # https://stackoverflow.com/questions/14848172/appending-a-list-to-a-list-of-lists-in-r
# # https://stackoverflow.com/questions/46310288/r-list-of-lists-to-dataframe-with-list-name-as-extra-column
#Info...
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6591304/
# https://stackoverflow.com/questions/55445629/appending-a-list-in-a-loop-r
# https://www.programmingr.com/fast-r-append-list/
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/nchar.html
# https://www.johnmyleswhite.com/notebook/2009/02/25/text-processing-in-r/
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/row.names.html
#Test Scina etc.
# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Automated_Cell_Type_Annotation
# https://satijalab.org/seurat/archive/v3.1/interaction_vignette.html
# https://datatofish.com/export-dataframe-to-csv-in-r/
##Enrichment test:
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
# https://bioconductor.org/packages/release/bioc/ html/fgsea.html
# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
#if you need more later:
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
