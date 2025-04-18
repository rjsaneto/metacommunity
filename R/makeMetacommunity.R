#' Create a Metacommunity in Near Neutral Communities Models
#'
#' Function to run the metacommunity to model in Near Neutral Metacommunities
#'
#'@param J The number of individuals in Species Pool
#'@param theta Hubbels Fundamental Number of Biodiversity
#'@param sigma the standard deviation for change the per capita fecundity parameter w
#'@param Jl The number of individuals in local species
#'@param nloc Number of local communities
#'@param mig Migration rate number
#'@param loccomm Option of get the per capita fecundity parameter in local communities: "pool", "uniform", "diversity" or "neutral"
#'@param abonly Logical Use only abundance of species pool in migration events instead include per capita fecundity parameter w
#'@param dataframe Logical output as data.frame and metacommunity object instead list object
#'
#'@returns metacomm function return a metacomm object or a data.frame with information of a set of local communities.
#'@returns Such Local communities could attributes a set of wi (e.g. loccom="diversity") or get a uniform wi for all local communities (e.g. loccom="diversity"), or keep the wi of species pool (e.g. loccom="pool").
#'@returns Also Local communities could be neutral if loccom="neutral" or if sigma = 0, in this options all metacommunity will be neutral.
#'@seealso iteration pool metacomm
#'@examples
#'#Keep the wi of specie pool
#'set.seed(2024)
#'myMeta<-metacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="pool",abonly=F,dataFrame=T)
#'myMeta #see the species pool
#'tapply(myMeta$wi,list(myMeta$species,myMeta$source),sum)
#'
#' # Uniform wi in Local community
#'set.seed(2024)
#'myMeta<-metacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="uniform",abonly=F,dataFrame=T)
#'myMeta #see the species pool
#'tapply(myMeta$wi,list(myMeta$species,myMeta$source),sum)
#'
#' # Diversity wi in Local community
#'set.seed(2024)
#'myMeta<-metacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="diversity",abonly=F,dataFrame=T)
#'myMeta #see the species pool
#'tapply(myMeta$wi,list(myMeta$species,myMeta$source),sum)
#'
#' # Neutral Local community
#'set.seed(2024)
#'myMeta<-metacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="neutral",abonly=F,dataFrame=T)
#'myMeta #see the species pool
#'tapply(myMeta$wi,list(myMeta$species,myMeta$source),sum)
#'
#' # Neutral metacommunity
#'set.seed(2024)
#'myMeta<-metacomm(J=1000,theta = 100,sigma=0,Jl=100,nloc=3,mig=0.25,loccomm="diversity",abonly=F,dataFrame=T)
#'myMeta #see the species pool
#'tapply(myMeta$wi,list(myMeta$species,myMeta$source),sum)
#'@export
metacomm<-function(J,theta,sigma,Jl,nloc,mig,loccomm="pool",abonly=F,dataFrame=T,exp=T,itera=T){
  # function to generate a metacommunity
  mypool<-pool(J = J,theta = theta,sigma = sigma,ord = T,exp=exp)
  if(itera){
    mypool<-iteration(mypool,breaks = 1e6,exp=exp)
  }

  mylocal<-list()
  ## jl of the same number of individual, but we can insert a random fission procedure
  jl<-round(Jl/nloc,digits = 0)
  jl<-rep(jl,nloc)

  sigloc=NULL
  wi=NULL
  if(loccomm=="uniform"){
    if(exp){
      wi=rnorm(n = nrow(mypool),mean = 0,sd = sigma)
    }else{
      wi=rnorm(n = nrow(mypool),mean = 1,sd = sigma)
    }
  }

  if(loccomm=="diversity"){
    sigloc=sigma
  }

  if(loccomm=="neutral"){
    sigloc=0
  }

  for(i in 1:length(jl)){
    mylocal[[i]]<-one.local(pool = mypool,jl = jl[i],mig = mig,sigma = sigloc,wi = wi,abonly = abonly,exp = exp)
    if(itera){
      mylocal[[i]]<-iteration(mylocal[[i]],breaks = 1e6,exp = exp)
    }
  }

  if(dataFrame){
    results<-as.data.frame(mypool)
    results$source<-"pool"
    for(i in 1:length(mylocal)){
      temp<-as.data.frame(mylocal[[i]])
      temp$source<-paste("local.",i,sep="")
      results<-rbind(results,temp)
    }
    results$sigma<-sigma
    results$feature<-loccomm
    results$theta<-theta
    if(abonly){
      results$migration<-"neutral"
    }else{
      results$migration<-"fitness"
    }
    class(results)<-c("metacomm","data.frame")
  }else{
    results<-mylocal
    names(results)<-paste("local.",1:length(results),sep="")
    results$pool<-mypool
  }
  return(results)
}

#'@export
print.metacomm<-function(x){
  temp<-cbind((tapply(x$abundance,INDEX = list(x$source),FUN = richness)),(tapply(x$abundance,INDEX = list(x$source),FUN = sum)))
  colnames(temp)<-c("richness","abundance")
  print(t(temp))
}
