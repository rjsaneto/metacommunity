#'@export
one.local<-function(pool,jl,mig,sigma=NULL,wi=NULL,abonly=F){
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
    comm$wi<-rnorm(n = nrow(comm),mean = 1,sd = sigma)
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
      probPool<-rad(pool$abundance,pool$wi)
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
