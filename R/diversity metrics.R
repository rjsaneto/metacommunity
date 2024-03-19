richness<-function(abundance){
  rich<-length(abundance[abundance>0])
  return(rich)
}

rad<-function(abundance,w=NULL){
  if(!is.null(w)){
    abundance<-abundance*w
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
