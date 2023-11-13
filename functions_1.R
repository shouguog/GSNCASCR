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
  metaData_type1_perm<-metaData_type1type2[1:length(metaData_type1[,1]),]
  metaData_type2_perm<-metaData_type1type2[1:length(metaData_type1[,1])+length(metaData_type2[,1]),]
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

# Evaluate significance of each set by comparing
# to set scores after permutation of group labels
network_seq_correlation = function(pathway,seuratObj,genes_selected,type1,type2, returnType="network"){
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
  CSCORE_type1_coexp<-CSCORE_type1_coexp[rownames(CSCORE_type1_coexp) %in% pathway, colnames(CSCORE_type1_coexp) %in% pathway]
  CSCORE_type2_coexp<-CSCORE_type2_coexp[rownames(CSCORE_type2_coexp) %in% pathway, colnames(CSCORE_type2_coexp) %in% pathway]
  list(CSCORE_type1_coexp=CSCORE_type1_coexp, CSCORE_type2_coexp=CSCORE_type2_coexp)
}

network_plot = function(networklist,name1="status1",name2="status2", fileName=NULL){
  CSCORE_type1_coexp<-networklist[[1]]
  CSCORE_type2_coexp<-networklist[[2]]
  ####Need copy edit
  library(ggraph)
  gene1<-c()
  gene2<-c()
  r1<-c()
  r2<-c()
  for(ii in 1:(dim(CSCORE_type1_coexp)[1]-1)){
    for(jj in (ii+1):dim(CSCORE_type1_coexp)[1]){
      gene1<-c(gene1, rownames(CSCORE_type1_coexp)[ii])
      gene2<-c(gene2, rownames(CSCORE_type1_coexp)[jj])
      r1<-c(r1, as.numeric(CSCORE_type1_coexp[ii,jj]))
      r2<-c(r2, as.numeric(CSCORE_type1_coexp[ii,jj]))
    }
  }
  network_cors_1<-data.frame(x=gene1, y=gene2, r1=r1, r2=r2)
  graph_cors <- network_cors_1 %>% filter((abs(r1)+abs(r2)) > .6) %>% graph_from_data_frame(directed = FALSE)
  
  #and here’s one that’s polished to look nicer:
  gg1<-ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r1), edge_width = abs(r2), color = r1)) + guides(edge_alpha = "none", edge_width = "none") 
  gg1<-gg1 + scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) + geom_node_point(color = "white", size = 5) 
  gg1<-gg1 + geom_node_text(aes(label = name), repel = TRUE) + labs(title = name1)
  
  gg2<-ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r2), edge_width = abs(r2), color = r2)) + guides(edge_alpha = "none", edge_width = "none") 
  gg2<-gg2 + scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) + geom_node_point(color = "white", size = 5) 
  gg2<-gg2 + geom_node_text(aes(label = name), repel = TRUE) + labs(title = name2) 
  
  ggDiff<-ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r2), edge_width = abs(r2), color = (r2-r1))) + guides(edge_alpha = "none", edge_width = "none") 
  ggDiff<-ggDiff + scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) + geom_node_point(color = "white", size = 5) 
  ggDiff<-ggDiff + geom_node_text(aes(label = name), repel = TRUE) + labs(title = "Difference") 
  if(!is.null(fileName)){
    png(paste0(fileName, ".png"), width = 9000, height = 3000, res = 200)
    gridExtra::grid.arrange(gg1, gg2, gg3, ncol=3)
    dev.off()
  }
  return(list(gg1=gg1, gg2=gg2, ggDiff=ggDiff))
  ####Need copy edit
}


####The main program of runPermutation
####ncores control the pararel or not
####Seems that parael can not work in windows sometime
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


