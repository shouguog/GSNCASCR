##R 4.1.2

library(CSCORE)
library(Seurat)
library(parallel)
library(doParallel)
library(GSEABase)
## compute set correlation #
setsCorrelation = function(setInd,corDiff){
  setInd<-setInd[setInd %in% rownames(corDiff)]
  #cat(setInd)
  corDiff[setInd, setInd]
}
## compute set correlation scores. Now used sum function
##In the future extend to GSNCA score
setsCorrelationScore = function(setInd,corDiff){
  setInd<-setInd[setInd %in% rownames(corDiff)]
  sum(abs(corDiff[setInd, setInd]))
}

## Calculate with two networks
GSNCAscore = function(cormat1, cormat2){
  cormat1 <- abs(cormat1) #abs(cor(group1, method=cor.method))
  cormat2 <- abs(cormat2) #abs(cor(group2, method=cor.method))
  e1 <- eigen(cormat1)
  e2 <- eigen(cormat2)
  p1 <- abs(e1$vectors[,1])
  p2 <- abs(e2$vectors[,1])
  D_obs <- sum(abs((p1*norm(matrix(p1))) - (p2*norm(matrix(p2)))))
  D_obs
}

# Evaluate significance of each set by comparing
# to set scores after permutation of group labels
doPermute_seq_correlation = function(sets,seuratObj,genes_selected,type1,type2, returnType="network", seedPerm=0){
  metaData<-seuratObj@meta.data
  metaData_type1<-metaData[metaData$Status==type1,]
  metaData_type2<-metaData[metaData$Status==type2,]
  ###Get the metaData
  if(length(metaData_type1)!=length(metaData_type1)){
    stop("Please keep the sample equal size\n")
  }
  metaData_type1type2<-rbind(metaData_type1, metaData_type2)
  ###Start permutation
  set.seed(seedPerm)
  metaData_type1type2<-metaData_type1type2[sample(1:dim(metaData_type1type2)[1], dim(metaData_type1type2)[1]),]
  metaData_type1_perm<-metaData_type1type2[1:length(metaData_type1),]
  metaData_type2_perm<-metaData_type1type2[1:length(metaData_type1)+length(metaData_type2),]
  cat(head(rownames(metaData_type1_perm))[1:5])
  cat(head(rownames(metaData_type2_perm))[1:5])
  ###Get Seurat with samples
  pbmc_B_type1 <- subset(seuratObj, cells = rownames(metaData_type1_perm))
  pbmc_B_type2 <- subset(seuratObj, cells = rownames(metaData_type2_perm))
  
  CSCORE_result_type1 <- CSCORE(pbmc_B_type1, genes = genes_selected)
  CSCORE_result_type2 <- CSCORE(pbmc_B_type2, genes = genes_selected)
  
  # Obtain CS-CORE co-expression estimates
  CSCORE_type1_coexp <- CSCORE_result_type1$est
  # Obtain CS-CORE co-expression estimates
  CSCORE_type2_coexp <- CSCORE_result_type2$est
  corDiff=CSCORE_type2_coexp-CSCORE_type1_coexp
  #cat("I am here line 62\n")
  if(returnType=="network"){
    sapply(sets,setsCorrelation,corDiff=corDiff)
  }else{
    sapply(sets,setsCorrelationScore,corDiff=corDiff)
  }
}

# Evaluate significance of each set by comparing
# to set scores after permutation of group labels
calculate_seq_correlation = function(sets,seuratObj,genes_selected,type1,type2, returnType="network"){
  ###Get the metaData
  if(length(type1)!=length(type2)){
    stop("Please keep the sample equal size\n")
  }
  metaData<-seuratObj@meta.data
  metaData_type1<-metaData[metaData$Status==type1,]
  metaData_type2<-metaData[metaData$Status==type2,]
  ###Get Seurat with samples
  pbmc_B_type1 <- subset(seuratObj, cells = rownames(metaData_type1))
  pbmc_B_type2 <- subset(seuratObj, cells = rownames(metaData_type2))
  
  CSCORE_result_type1 <- CSCORE(pbmc_B_type1, genes = genes_selected)
  CSCORE_result_type2 <- CSCORE(pbmc_B_type2, genes = genes_selected)
  
  # Obtain CS-CORE co-expression estimates
  CSCORE_type1_coexp <- CSCORE_result_type1$est
  # Obtain CS-CORE co-expression estimates
  CSCORE_type2_coexp <- CSCORE_result_type2$est
  corDiff=CSCORE_type2_coexp-CSCORE_type1_coexp
  cat("I am here line 62\n")
  if(returnType=="network"){
    sapply(sets,setsCorrelation,corDiff=corDiff)
  }else{
    sapply(sets,setsCorrelationScore,corDiff=corDiff)
  }
}

####The main program of runPermutation
####ncores control the pararel or not
####Seems that parael can not work in windows sometime, need to fix in future
####Need to redefine function in when ncores>1, need to find other solution
runPermute_seq_correlation = function(sets,seuratObj,genes_selected,type1,type2, ncores=1,nPermute=6, returnType="network"){
  if(ncores==1){
    correlation_perm_list<-list()
    for(ii in 1:nPermute){
      cat("############################################");cat(Sys.time());cat("\t");cat(ii); cat("\tperm start\n")
      correlation_perm_list[[ii]]<-doPermute_seq_correlation(sets,seuratObj,genes_selected,type1,type2, returnType=returnType, seedPerm=ii)
    }
    correlation_perm_list
  }else{
    cl = makeCluster(ncores)
    # register the parallel backend
    registerDoParallel(cl)
    # Set path to personal library on all child processes
    clusterEvalQ(cl,.libPaths('~/Rlib'))
    cat('Starting run 1 ...\n')
    cat(dim(seuratObj@meta.data)[1]);cat("\n")
    correlation_perm_list = foreach(n=1:nPermute) %dopar% {
      ##In the future extend to GSNCA score
    setsCorrelationScore = function(setInd,corDiff){
      setInd<-setInd[setInd %in% rownames(corDiff)]
      sum(abs(corDiff[setInd, setInd]))
    }
    
    # Evaluate significance of each set by comparing
    # to set scores after permutation of group labels
    doPermute_seq_correlation = function(sets,seuratObj,genes_selected,type1,type2, returnType="network", seedPerm=0){
      library(Seurat)
      library(CSCORE)
      metaData<-seuratObj@meta.data
      metaData_type1<-metaData[metaData$Status==type1,]
      metaData_type2<-metaData[metaData$Status==type2,]
      ###Get the metaData
      if(length(metaData_type1)!=length(metaData_type1)){
        stop("Please keep the sample equal size\n")
      }
      metaData_type1type2<-rbind(metaData_type1, metaData_type2)
      ###Start permutation
      set.seed(NULL)
      metaData_type1type2<-metaData_type1type2[sample(1:dim(metaData_type1type2)[1], dim(metaData_type1type2)[1]),]
      metaData_type1_perm<-metaData_type1type2[1:length(metaData_type1),]
      metaData_type2_perm<-metaData_type1type2[1:length(metaData_type1)+length(metaData_type2),]
      cat(head(rownames(metaData_type1_perm))[1:5])
      cat(head(rownames(metaData_type2_perm))[1:5])
      ###Get Seurat with samples
      pbmc_B_type1 <- subset(seuratObj, cells = rownames(metaData_type1_perm))
      pbmc_B_type2 <- subset(seuratObj, cells = rownames(metaData_type2_perm))
      
      CSCORE_result_type1 <- CSCORE(pbmc_B_type1, genes = genes_selected)
      CSCORE_result_type2 <- CSCORE(pbmc_B_type2, genes = genes_selected)
      
      # Obtain CS-CORE co-expression estimates
      CSCORE_type1_coexp <- CSCORE_result_type1$est
      # Obtain CS-CORE co-expression estimates
      CSCORE_type2_coexp <- CSCORE_result_type2$est
      corDiff=CSCORE_type2_coexp-CSCORE_type1_coexp
      #cat("I am here line 62\n")
      if(returnType=="network"){
        sapply(sets,setsCorrelation,corDiff=corDiff)
      }else{
        sapply(sets,setsCorrelationScore,corDiff=corDiff)
      }
    }
    cat(dim(seuratObj@meta.data)[1]);cat("\n")
    doPermute_seq_correlation(sets,seuratObj,genes_selected,type1,type2, returnType=returnType)
    }
    correlation_perm_list
  }
}

####Get the pvalues
getPvalues = function(cor_data_score,cor_data_perm_list){
  pvalues<-c()
  for(pathway in names(cor_data_score)){
    score<-cor_data_score[pathway]
    scores<-c()
    for(ii in 1:length(cor_data_perm_list)){
      scores<-c(scores, cor_data_perm_list[[ii]][pathway])
    }
    pvalues<-c(pvalues, sum(score>scores)/length(cor_data_perm_list))
  }
  names(pvalues)<-names(cor_data_score)
  pvalues
}


