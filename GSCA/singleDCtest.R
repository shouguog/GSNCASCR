rm(list = ls())

singleDC <- function(data, group, GSdefList, nperm){
  
  DI=rep(0,length(GSdefList))
  pvalue=rep(0,length(GSdefList))
  permv=matrix(0,length(GSdefList),nperm)
  
  for(j in 1:length(GSdefList)){
    permv[j,1]=compute.test.stat(data[GSdefList[[j]],],group)
  }
  
  print('Perm 1 completed')
  
  for(i in 2:nperm){
    data.perm=data[,sample(length(data[1,]))]
    
    for(j in 1:length(GSdefList)){
      permv[j,i]=compute.test.stat(data.perm[GSdefList[[j]],],group)
    }
    
    print(paste('Perm',i,'completed'))
  }
  
  DI=permv[,1]
  pvalue=rowSums(permv>=permv[,1])/length(permv[1,])
  
  ret=NULL
  
  names(DI)=names(GSdefList)
  names(pvalue)=names(GSdefList)
  dimnames(permv)[[1]]=names(GSdefList)
  permnames=NULL
  for(m in 1:nperm){
    permnames[m]=paste('P',m,sep='')
  }
  dimnames(permv)[[2]]=permnames
  
  ret$DI=DI
  ret$pvalue=pvalue
  ret$permv=permv
  
  return(ret)
}


compute.test.stat = function(fixed.gs.data,group){
  
  n.pairs<-choose(nrow(fixed.gs.data),2)
  n.groups<-length(group)
  csum.group<-cumsum(group)
  
  corr.mats=NULL
  
  if(n.groups==2){
    corr.mats[[1]]=cor(t(fixed.gs.data[,1:(group[1])]))
    corr.mats[[2]]=cor(t(fixed.gs.data[,(group[1]+1):(csum.group[2])]))}
  
  if(n.groups==3){
    corr.mats[[1]]=cor(t(fixed.gs.data[,1:(group[1])]))
    corr.mats[[2]]=cor(t(fixed.gs.data[,(group[1]+1):(csum.group[2])]))
    corr.mats[[3]]=cor(t(fixed.gs.data[,(csum.group[2]+1):(csum.group[3])])) }
  
  if(n.groups>3){
    for(i in 1:n.groups){
      if(i==1){corr.mats[[i]]=cor(t(fixed.gs.data[,1:group[i]]))} 
      else{corr.mats[[i]]=cor(t(fixed.gs.data[,(csum.group[i-1]+1):csum.group[i]]))}
    }}
  
  vals=NULL
  k=1
  for(i in 1:(n.groups-1)){
    for(j in 2:n.groups){
      sq.diff.mat=(corr.mats[[i]]-corr.mats[[j]])^2
      vals[[k]]=sqrt((1/n.pairs)*sum(sq.diff.mat[upper.tri(sq.diff.mat)]))
      k=k+1
    }
  }
  
  test.stat=mean(unlist(vals))  #shouguo test.stat=mean(vals)
  return(test.stat)
}               



load("data/LungCancer3.rda")
str(LungCancer3$data)
GS <- LungCancer3$info$GSdef
GS[["GO:0008033"]]

dc.M <- singleDC(data = LungCancer3$data$Michigan, group = c(86,10), GSdefList = GS, nperm = 3)
