##R 4.1.2
rm(list = ls())
library(CSCORE)
library(Seurat)
library(parallel)
library(doParallel)
library(GSEABase)
source("functions_0.R")
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

####Calculate the network, or score of networks
#cor_data<-calculate_seq_correlation(geneset,pbmc_B.sub,genes_selected,"COVID","Healthy")
cor_data_score<-calculate_seq_correlation(geneset,pbmc_B.sub,genes_selected,"COVID","Healthy", returnType = "score")
save(cor_data_score, file = "cor_data_score.RData")

####Calculate the network, or score of networks for permutations
cor_data_perm_list<-runPermute_seq_correlation(geneset,pbmc_B.sub,genes_selected,"COVID","Healthy", ncores = 1, nPermute=1000, returnType = "score")
save(cor_data_perm_list, file = "cor_data_hallmark_perm_list_1000.RData")


#cor_data_perm_list<-runPermute_seq_correlation(geneset,pbmc_B.sub,genes_selected,"COVID","Healthy", ncores = 4, nPermute=6, returnType = "score")

getPvalues(cor_data_score, cor_data_perm_list)

