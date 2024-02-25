#'@export
pool<-function(J,theta,sigma,ord=T){
  #initial attributes of species pool
  abundance<-1
  species<-1
  w<-rnorm(1,1,sigma)

  #pool creation procedure
  for(j in 2:J){
    dice<-runif(1)
    if(dice<(theta/(theta+j-1))){
      #if dice < theta... the new individual belong of a new species
      species<-c(species,max(species)+1)
      abundance<-c(abundance,1)
      w<-c(w,rnorm(1,1,sigma))
    }else{
      #else the new individual belong of one of already species of pool
      probtemp<-rad(abundance = abundance,w = w)

      who<-sample(x = 1:length(species),size = 1,prob = probtemp)
      abundance[who]<-abundance[who]+1

    }
  }

  #output production
  species<-paste("sp.",species,sep="")
  results<-as.data.frame(species)
  results$abundance<-abundance
  results$wi<-w
  if(ord){
    results<-results[order(abundance,decreasing = T),]
  }
  attr(results,"richness")<-richness(abundance)
  attr(results,"extinction")<-0
  attr(results,"theta")<-theta
  attr(results,"sigma")<-sigma
  class(results)<-c("pool","data.frame")
  return(results)
}

#'@export
print.pool<-function(x){
  cat("number of species:",richness(x$abundance),"\n")
  if(nrow(x)>10){
    print(head(as.data.frame(x),3))
    cat(rep("\t.\n",3))
    print(tail(as.data.frame(x),3))
  }else{
    print.data.frame(x)
  }
  invisible(x)
}
