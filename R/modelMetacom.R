#' Model a Metacommunity in Near Neutral Communities
#'
#' Function to run the model of Near Neutral Metacommunities
#'
#'@param J The number of individuals in Species Pool
#'@param theta Hubbels Fundamental Number of Biodiversity
#'@param sigma the standard deviation for change the per capita fecundity parameter w
#'@param Jl The number of individuals in local species
#'@param nloc Number of local communities
#'@param mig Migration rate number
#'@param loccomm Option of get the per capita fecundity parameter in local communities: "pool", "uniform", "diversity" or "neutral"
#'@param abonly Logical Use only abundance of species pool in migration events instead include per capita fecundity parameter w
#'@param nSim Number of simulations of the model
#'@param file Optional file to exporte the output to *.rData file
#'
#'@returns The model of metacommunity
#'@seealso iteration pool metacomm one.local
#'@examples
#'#Get an object in R
#'set.seed(2024)
#'
#'myModel<-modelMetacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="diversity",abonly=F,nSim=2)
#'myModel #see the species pool
#'
#' #Save in external file
#'set.seed(2024)
#'modelMetacomm(J=1000,theta = 100,sigma=0.03,Jl=100,nloc=3,mig=0.25,loccomm="diversity",abonly=F,nSim=2,file="mytest.RData")
#'load("mytest.RData") #see the species pool
#'model_1000_100_0.03_100_3_0.25_diversity

#'@export
modelMetacomm<-function(J,theta,sigma,Jl,nloc,mig,loccomm="pool",abonly=F,nSim=5,file=NULL){
  results<-NULL
  for(i in 1:nSim){
    temp<-metacomm(J=J,theta=theta,sigma=sigma,Jl=Jl,nloc=nloc,mig=mig,loccomm=loccomm,abonly=abonly,dataFrame=T)
    temp$simulation<-i
    results<-rbind(results,temp)
  }
  if(!is.null(file)){
    nome<-paste("model_",as.integer(J),"_",theta,"_",sigma,"_",Jl,"_",nloc,"_",mig,"_",loccomm,sep="")
    assign(nome,value = results)
    save(list = nome,file = file)
  }
  class(results)<-c("metacomm","data.frame")
  invisible(results)
}

#'@export
print.metacomm<-function(x){
  temp<-tapply(x$abundance, list(x$simulation,x$source), richness)
  media<-round(apply(X = temp,FUN = mean,MARGIN = 2),2)
  desvio<-round(apply(X = temp,FUN = sd,MARGIN = 2),2)
  cat("\n Species Richness of",max(x$simulation),"simulations\n\n")
  cat("source","\t\t--\t","mean","\t+/-  ","standard deviation","\n")
  for(i in 1:(length(media)-1)){
    cat(names(media)[i],"\t--\t ",media[i],"\t+/-  ",desvio[i],"\n")
  }
  cat(names(media)[length(media)],"\t\t--\t",media[length(media)],"\t+/- ",desvio[length(media)],"\n")
}
