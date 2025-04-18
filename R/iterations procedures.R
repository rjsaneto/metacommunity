#one.iteration class

one.iteration<-function(ob,exp=T){
  UseMethod("one.iteration")
}

one.iteration.default<-function(ob,exp=T){
  cat("The object is not a metacommunity")
}

#one.iteration for class pool
one.iteration.pool<-function(comm,exp=T){
  #this function will create a one iteration procedure of death-birth zero-sum
  #comm: community make of pool class

  rich<-nrow(comm)
  sp<-rich+1

  theta<-attr(comm,"theta")
  sigma<-attr(comm,"sigma")

  #death procedure
  probtemp<-rad(comm$abundance)
  who<-sample(x = 1:nrow(comm),size = 1,prob = probtemp)
  comm$abundance[who]<-comm$abundance[who]-1
  #birth procedure
  dice<-runif(1)
  ab<-sum(comm$abundance)
  if(dice<(theta/(theta+ab-1))){ #I suppose theta+ab and not theta+ab-1 once the remotion of individual already happens in death procedure
  #chance of new species
    #new row for new species
    comm<-rbind(comm,rep(NA,ncol(comm)))
    #fill the species column
    comm$species[sp]<-paste("sp.",sp,sep = "")
    #fill abundance column
    comm$abundance[sp]<-1
    #fill wi column
    if(exp){
      comm$wi[sp]<-rnorm(n = 1,mean = 0,sd = sigma)
    }else{
      comm$wi[sp]<-rnorm(n = 1,mean = 1,sd = sigma)
    }
    #increase +1 sp
    sp<-sp+1
  }else{
  #else another species of comm
    probtemp<-rad(comm$abundance,comm$wi,exp=exp)
    who<-sample(x = 1:nrow(comm),size = 1,prob = probtemp)
    comm$abundance[who]<-comm$abundance[who]+1
  }
  attr(comm,"extinction")<-c(attr(comm,"extinction"),length(comm$abundance[comm$abundance==0]))
  attr(comm,"richness")<-c(attr(comm,"richness"),richness(comm$abundance))
  return(comm)
}

#one.iteration for class local
#'@export
one.iteration.local<-function(comm,exp=T){
  #this function will create a one iteration procedure of death-birth zero-sum in a local community
  #comm: local community make of local class

  #extract attributes
  rich<-nrow(comm)

  probpool<-attr(comm,"probpool")
  sigma<-attr(comm,"sigma")
  mig<-attr(comm,"migration")

  #death procedure
  probtemp<-rad(comm$abundance)
  who<-sample(x = 1:nrow(comm),size = 1,prob = probtemp)
  comm$abundance[who]<-comm$abundance[who]-1

  #bird procedure
  dice<-runif(1)
  ab<-sum(comm$abundance)
  if(dice<mig){ #I could implement (mig/(mig+ab) and not mig but later
    #chance of species from pool
    who<-sample(x=1:nrow(comm),size=1,prob = probpool)
    comm$abundance[who]<-comm$abundance[who]+1
  }else{
    #else another species of comm
    probtemp<-rad(comm$abundance,comm$wi,exp=exp)
    who<-sample(x = 1:nrow(comm),size = 1,prob = probtemp)
    comm$abundance[who]<-comm$abundance[who]+1
  }
  ex<-attr(comm,"local.extinction")
  ex<-ex[length(ex)]

  rich<-richness(comm$abundance)
  richold<-attr(comm,"richness")
  richold<-richold[length(richold)]
  if(rich<richold){
    ex<-ex+1
  }
  attr(comm,"local.extinction")<-c(attr(comm,"local.extinction"),ex)
  attr(comm,"richness")<-c(attr(comm,"richness"),rich)
  return(comm)
}

#######################################################################################
#' Iteration procedure of Near Neutral Communities
#'
#' Function to run a iteration procedure in Local and Pool Communities
#'
#'@param comm Object of class pool or local
#'@param breaks Maximum number of iteration
#'
#'@returns Iteration function return the original class object after run several iterations of death and birth events, such events should be accessed by attributes,see example.
#'@seealso pool one.local metacomm
#'@examples
#'set.seed(2024)
#'mypool<-pool(J=1000,theta = 100,sigma=0.03)
#'mypool #see the species pool
#'attr(mypool,"richness") #only one richness recorded
#'
#'mynewpool<-iteration(mypool)
#'mynewpool
#'attr(mynewpool,"richness") #a set of richness recorded
#'@export
iteration<-function(comm,breaks=1e6,exp=T){
  i=0

  #isso vai sair
  medRiq<-NULL
  varRiq<-NULL

  #matriz
  commMatrix<-extractComm(comm)
  commMatrix<-commMatrix/sum(commMatrix)

  repeat{

    comm<-one.iteration(comm,exp=exp)
    new<-extractComm(comm)
    new<-new/sum(new)

    commMatrix<-combineMatrix(new,commMatrix)

    riq<-attr(comm,"richness")
    varRiq<-c(varRiq,mean(riq))

    i<-i+1
    if(i>2000){
      medRiq<-c(medRiq,mean(riq[(i-2000):i]))

      med1<-mean(riq)
      med2<-mean(riq[1:(length(riq)-1)])

      var1<-var(riq)
      var2<-var(riq[1:(length(riq)-1)])

      diferMed<-abs(med1-med2)
      diferVar<-abs(sqrt(var1)-sqrt(var2))
      diferRate<-abs((var1/med1)-(var2/med2))

      if(breaks<=length(riq)){
        break
      }
      if(abs(calcDiff(commMatrix))<=1e-18){

        if(abs(medRiq[length(medRiq)]-riq[length(riq)])<1){
          cat(paste(rep("#",20)),paste("\n ####     ENCERRADO!!!! ##### \n"),paste(rep("#",20)),paste("\n"))
          break
        }else{
          cat(
            paste("Diferença =",abs(calcDiff(commMatrix)),"menor que requerida para a Iteração:", i," Riqueza Media =",medRiq[length(medRiq)]," Riqueza Atual = ",riq[length(riq)],"\n"))
        }
      }else{
        cat(paste("Diferença =",abs(calcDiff(commMatrix)),"maior que requerida para a Iteração:", i,"\n"))
      }

    }else{cat(paste("Iniciando Iteração:", i,"\n"))}

  }
  attr(comm,"variance")<-varRiq
  attr(comm,"average")<-medRiq
  attr(comm,"matriz")<-commMatrix
  return(comm)
}
