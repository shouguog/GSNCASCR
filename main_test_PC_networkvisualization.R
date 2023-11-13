##R 4.1.2
rm(list = ls())
library(CSCORE)
library(Seurat)
library(parallel)
library(doParallel)
library(GSEABase)
source("functions_1.R")
####Load the gene set and convert into list
geneset<-geneIds(getGmt("h.all.v2023.2.Hs.symbols.gmt"))
####Load the Seurat
load("pbmc_B.RData")
metaData<-pbmc_B@meta.data
####make sure two groups have same number of cells for comparisiion
metaData_COVID<-metaData[metaData$Status=="COVID",]
metaData_Healthy<-metaData[metaData$Status=="Healthy",]
set.seed(666);metaData_COVID<-metaData_COVID[sample(1:dim(metaData_COVID)[1],dim(metaData_Healthy)[1]),]
pbmc_B.sub<-subset(pbmc_B, cells=c(rownames(metaData_COVID), rownames(metaData_Healthy)))
rm(pbmc_B)

for(pathway in names(geneset)[c(1,25,27,45)]){
  networkslist<-network_seq_correlation(geneset[[pathway]],pbmc_B.sub,genes_selected,"COVID","Healthy")
  gglist<-network_plot(networkslist,name1="COVID",name2="Healthy", fileName=paste0(pathway, "_COVID_Healthy_Diff"))
}
