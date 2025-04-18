richness<-function(abundance){
  rich<-length(abundance[abundance>0])
  return(rich)
}

rad<-function(abundance,w=NULL,exp=T){
  if(!is.null(w)){
    if(exp){
      abundance<-abundance*exp(w)
    }else{
      abundance<-abundance*w
    }
  }
  abundance<-abundance/sum(abundance)
  return(abundance)
}

cummean<-function(x){
  result<-NULL
  for(i in 1:length(x)){
    result<-c(result,mean(x[1:i]))
  }
  return(result)
}

extractComm<-function(df){
  results<-df$abundance
  results<-results[results>0]
  results<-results[order(results,decreasing = T)]
}

combineMatrix<-function(new,old){
  if(!is.matrix(old)){
    old<-as.matrix(old)
  }
  if(length(new)>nrow(old)){
    zeros<-matrix(0,nrow = length(new)-nrow(old),ncol = ncol(old))
    old<-rbind(old,zeros)
  }
  old<-cbind(old,0)
  old[1:length(new),ncol(old)]<-new
  colnames(old)<-paste("iter.",1:ncol(old),sep="")
  return(old)
}

calcDiff<-function(mat){
  newMean<-rowMeans(mat)
  oldMean<-rowMeans(mat[,-ncol(mat)])
  differ<-newMean-oldMean
  return(mean(differ))
}
