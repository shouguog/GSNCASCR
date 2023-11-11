# This is need to be modified. 
# https://stackoverflow.com/questions/20704235/function-not-found-in-r-doparallel-foreach-error-in-task-1-failed-cou
# Evaluate significance of each set by comparing
# to set scores after permutation of group labels
calculate_seq_correlation_para = function(sets,seuratObj,genes_selected,type1,type2, returnType="network"){
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

# to set scores after permutation of group labels
doPermute_seq_correlation_para = function(sets,seuratObj,genes_selected,type1,type2, returnType="network"){
  metaData<-seuratObj@meta.data
  metaData_type1<-metaData[metaData$Status==type1,]
  metaData_type2<-metaData[metaData$Status==type2,]
  ###Get the metaData
  if(length(metaData_type1)!=length(metaData_type1)){
    stop("Please keep the sample equal size\n")
  }
  metaData_type1type2<-rbind(metaData_type1, metaData_type2)
  ###Start permutation
  metaData_type1type2<-metaData_type1type2[sample(1:dim(metaData_type1type2)[1], dim(metaData_type1type2)[1]),]
  metaData_type1_perm<-metaData_type1type2[1:length(metaData_type1),]
  metaData_type2_perm<-metaData_type1type2[1:length(metaData_type1)+length(metaData_type2),]
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

