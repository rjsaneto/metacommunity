#'@export
modelMetacomm<-function(J,theta,sigma,Jl,nloc,mig,loccomm="pool",abonly=F,nSim=5,file=NULL){
  results<-NULL
  for(i in 1:nSim){
    temp<-metacomm(J=J,theta=theta,sigma=sigma,Jl=Jl,nloc=nloc,mig=mig,loccomm=loccomm,abonly=abonly,dataFrame=T)
    temp$simulation<-i
    results<-rbind(results,temp)
  }
  if(!is.null(file)){
    nome<-paste("model_",J,"_",theta,"_",sigma,"_",Jl,"_",nloc,"_",mig,"_",loccomm,sep="")
    assign(nome,value = results)
    save(list = nome,file = file)
  }
  class(results)<-c("metacomm","data.frame")
  return(results)
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
