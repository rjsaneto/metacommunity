#' Create a local community in Near Neutral Communities Models
#'
#' Function to run the local community of species to model in Near Neutral Metacommunities
#'
#'@param pool An object of class pool
#'@param jl The number of individuals in local species
#'@param mig Migration rate number
#'@param sigma the standard deviation for change the per capita fecundity parameter w
#'@param wi Vector of wi to impose uniformity in wi for a set of local communities
#'@param abonly Logical Use only abundance of species pool in migration events instead include per capita fecundity parameter w
#'
#'@returns pool function return a local class object with Species label, Abundance and per capita fecundity parameter w. Such object also record other parameters in attributes
#'@seealso iteration pool metacomm
#'@examples
#'set.seed(2024)
#'mypool<-pool(J=1000,theta = 100,sigma=0.03)
#'mypool #see the species pool
#'attributes(mypool)
#'attr(mypool,"richness") #only one richness recorded
#'
#'mylocal<-one.local(mypool,jl=100,mig=0.3,sigma=0.02,wi=NULL,abonly=F)
#'attributes(mylocal)
#'attr(mylocal,"richness")

#'@export
one.local<-function(pool,jl,mig,sigma=NULL,wi=NULL,abonly=F,exp=T){
  #function to creation of local community
  comm<-as.data.frame(pool)
  comm$abundance<-0
  if(is.vector(wi)){
    if(length(wi)!=nrow(comm)){
      wi<-NULL
      warning("Wi vector dont have the same length of pool community, so this can't be exported")
    }
  }
  if(!is.null(sigma)&&!is.null(wi)){
    sigma<-NULL
    warning("Both Sigma and Wi parameters are not null, so export of local wi vector is used, to create own local community wi, wi parameters need to be null")
  }

  #### If sigma is not null, create local own wi
  if(!is.null(sigma)){
    if(exp){
      comm$wi<-rnorm(n = nrow(comm),mean = 0,sd = sigma)
    }else{
      comm$wi<-rnorm(n = nrow(comm),mean = 1,sd = sigma)
    }
  }

  ### If wi is not null and have the same length of comm, export the wi
  if(!is.null(wi)){
    comm$wi<-wi
  }

  ### If both parameters sigma and wi are null, use the same wi of pool and don't change original wi.

  #### if the there is no migration fitness use only abundance to sample probability
  if(abonly==T){
      probPool<-rad(pool$abundance)
  }else{
      probPool<-rad(pool$abundance,pool$wi,exp=exp)
  }

  ##### first individual from pool
  who<-sample(1:nrow(comm),size = 1,prob = probPool)
  comm$abundance[who]<-comm$abundance[who]+1

  #### sequence of local community
  for(i in 2:jl){
    dice<-runif(1)
    # put here the migration probability condition
    if(dice<mig){
      # if dice < mig the new species was from pool
      who<-sample(1:nrow(comm),size = 1,prob = probPool)
      comm$abundance[who]<-comm$abundance[who]+1
    }else{
      #else the new species was from another local species
      problocal<-rad(comm$abundance,comm$wi)
      who<-sample(1:nrow(comm),size = 1,prob = problocal)
      comm$abundance[who]<-comm$abundance[who]+1
    }
  }

  #### output production

  attr(comm,"richness")<-richness(comm$abundance)
  attr(comm,"local.extinction")<-0
  attr(comm,"migration")<-mig
  attr(comm,"sigma")<-sigma
  attr(comm,"pool")<-pool$abundance
  attr(comm,"probpool")<-probPool
  class(comm)<-c("local","data.frame")
  return(comm)
}

#'@export
print.local<-function(x,show.all=F){
  temp<-x[x$abundance>0,]
  cat("number of species:",richness(temp$abundance),"\n")
  if(show.all==T){
    print.data.frame(temp)
  }else{
    if(nrow(temp)>10){
      print(head(as.data.frame(temp),3))
      cat(rep("\t.\n",3))
      print(tail(as.data.frame(temp),3))
    }else{
      print.data.frame(temp)
    }
  }
}
