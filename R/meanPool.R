#' Calculate the relative species abundance and W of the pool.
#'
#' Function to extract the arithmetic mean of the RSA and W from all simulations performed by the modelMetacomm function.
#'
#' @param df A dataframe containing the results of the simulations.
#' @returns The meanPool function returns an object containing the species rank, RSA (Relative Species Abundance), RSA standard deviation, W, W standard deviation and the simulation's sigma value.
#' @examples
#' mymodel <- modelMetacomm(J=100,theta = 10,sigma=0.03,Jl=10,nloc=3,mig=0.25,loccomm="diversity",abonly=F,nSim=2)
#' results <- meanPool(mymodel)


#'@export
meanPool <- function(df) {
  library(dplyr)
  simPool <- df %>% filter(source=="pool")
  simPool <- select(simPool, simulation, abundance, wi)

  #agora ordenaremos na seguinte ordem: simulação (1 a x), a abundancia da simulação x será ordenada em decrescente e o wi acompanhará a posição especifica da sua abundância.
  simPool <- simPool[with(simPool, order(simulation, -abundance)),]

  #criaremos a matriz para receber as simulações
  matrizRSA <- matrix(0, nrow = max(table(simPool$simulation)),
                      ncol = max(simPool$simulation))

  #for loop para inserir as abundâncias de cada simulação do pool em uma coluna
  for (i in 1:ncol(matrizRSA)) {
    matrizRSA[1:length(simPool$abundance[simPool$simulation == i]), i] <- simPool$abundance[simPool$simulation == i]
  }

  #criaremos a mtriz para receber os valores de wi
  matrizWi <- matrix(0, nrow = max(table(simPool$simulation)),
                     ncol = max(simPool$simulation))

  #for loop para colocar os wi's de cada simulação em uma coluna
  for (z in 1:ncol(matrizWi)) {
    matrizWi[1:length(simPool$wi[simPool$simulation == z]), z] <- simPool$wi[simPool$simulation == z]

  }

  #colocando NA nos zeros
  matrizRSA[matrizRSA == 0] <- NA
  matrizWi[matrizWi == 0] <- NA
  #tirar e médias e desvios das abundâncias e do wi e o rank
  ABUNDANCIA <- rowMeans(matrizRSA, na.rm = T)
  SD_ABUNDANCIA <- apply(matrizRSA, 1, sd, na.rm=T)
  WI <- rowMeans(matrizWi, na.rm = T)
  SD_WI <- apply(matrizWi, 1, sd, na.rm=T)
  RANK <- 1:length(ABUNDANCIA)
  final <- cbind(RANK, ABUNDANCIA, SD_ABUNDANCIA, WI, SD_WI)
  final <- as.data.frame(final) #coloquei como dataframe pois ocorreu erro no próximo passo

  final <- final %>% filter(ABUNDANCIA != 0) #aqui estou retirando linhas cuja abundancia é 0
  #adicionar os valroes de sigma e migração
  final$SIGMA <- rep(df$sigma[1], length(final$ABUNDANCIA))
  final$MIG <- rep(df$migration[1], length(final$ABUNDANCIA))
  return(final)
}




