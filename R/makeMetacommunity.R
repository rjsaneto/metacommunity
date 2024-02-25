#'@export
metacomm<-function(J,theta,sigma,Jl,nloc,mig,loccomm="pool",abonly=F,dataFrame=T){
  # function to generate a metacommunity
  mypool<-pool(J = J,theta = theta,sigma = sigma,ord = T)
  mypool<-iteration(mypool,breaks = 1e6)

  mylocal<-list()
  ## jl of the same number of individual, but we can insert a random fission procedure
  jl<-round(Jl/nloc,digits = 0)
  jl<-rep(jl,nloc)

  sigloc=NULL
  if(loccomm=="uniform"){
    wi=rnorm(n = nrow(mypool),mean = 1,sd = sigma)
  }

  if(loccomm=="diversity"){
    wi=NULL
    sigloc=sigma
  }

  for(i in 1:length(jl)){
    mylocal[[i]]<-one.local(pool = mypool,jl = jl[i],mig = mig,sigma = sigloc,wi = wi,abonly = abonly)
    mylocal[[i]]<-iteration(mylocal[[i]],breaks = 1e6)
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
