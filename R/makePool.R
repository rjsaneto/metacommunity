#' Create a Pool of Species in Near Neutral Communities Models
#'
#' Function to run the pool of species to model in Near Neutral Metacommunities
#'
#'@param J The number of individuals in Species Pool
#'@param theta Hubbels Fundamental Number of Biodiversity
#'@param sigma the standard deviation for change the per capita fecundity parameter w
#'@param ord Logical to make the decreasing ordination of abundance in species pool
#'
#'@returns pool function return a pool class object with Species label, Abundance and per capita fecundity parameter w. Such object also record richness and initial parameters in attributes
#'@seealso iteration one.local metacomm
#'@examples
#'set.seed(2024)
#'mypool<-pool(J=1000,theta = 100,sigma=0.03)
#'mypool #see the species pool
#'attributes(mypool)
#'attr(mypool,"richness") #only one richness recorded
#'
#'mynewpool<-iteration(mypool)
#'mynewpool
#'attr(mynewpool,"richness") #a set of richness recorded

#'@export
pool<-function(J,theta,sigma,ord=T,exp=T){
  #initial attributes of species pool
  abundance<-1
  species<-1
  if(exp){
    w<-rnorm(1,0,sigma)
  }else{
    w<-rnorm(1,1,sigma)
  }

  #pool creation procedure
  for(j in 2:J){
    dice<-runif(1)
    if(dice<(theta/(theta+j-1))){
      #if dice < theta... the new individual belong of a new species
      species<-c(species,max(species)+1)
      abundance<-c(abundance,1)
      if(exp){
        w<-c(w,rnorm(1,0,sigma))
      }else{
        w<-c(w,rnorm(1,1,sigma))
      }
    }else{
      #else the new individual belong of one of already species of pool
      #probtemp<-rad(abundance = abundance,w = rep(1,length(abundance)))
      probtemp<-rad(abundance = abundance,w = w,exp=exp)

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
  temp<-x[x$abundance>0,]
  temp<-temp[order(temp$abundance,decreasing=T),]
  if(nrow(temp)>10){
    print(head(as.data.frame(temp),3))
    cat(rep("\t.\n",3))
    print(tail(as.data.frame(temp),3))
  }else{
    print.data.frame(temp)
  }
  invisible(x)
}
