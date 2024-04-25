#' A function to calculate the relative species abundance of the local communities.
#'
#' Function to extract the arithmetic mean of the RSA from all simulations performed by the modelMetacomm function.
#'
#' @param df A dataframe containing the results of the simulations.
#' @param meanRSA If set to TRUE, the RSA of the simulations will be used for the mean; if set to FALSE, the RSA will be calculated from the mean abundances of all simulations. Generally, the first method may result in values different from 1 when the RSA is summed, whereas in the second method, this does not occur.
#' @returns The meanLocal function returns an object containing the species rank, RSA (Relative Species Abundance), RSA standard deviation, and the simulation's sigma value.
#' @examples
#' mymodel <- modelMetacomm(J=100,theta = 10,sigma=0.03,Jl=10,nloc=3,mig=0.25,loccomm="diversity",abonly=F,nSim=2)
#' results <- meanLocal(mymodel, meanRSA = TRUE) #The final RSA will be obtained by the means of the RSAs of all simulations.
#' results2 <- meanLocal(mymodel, meanRSA = FALSE) #The final RSA will be obtained by the  mean of the abundances from all simulations.


#'@export
meanLocal <- function(df, meanRSA = TRUE) {
  library(dplyr)
  #Filtrando o data frame apenas para local
  simLocais <- df %>% filter(source=="local.1")

  #criando matriz para abrigar as simulações realizadas
  matrizLocais <- matrix(0, nrow = max(table(simLocais$simulation)), #pega o número máximo de linhas da simulação com mais linhas
                         ncol = max(simLocais$simulation)) #pega o número máximo de simulações

  if (meanRSA == TRUE) {
    for (i in 1: max(simLocais$simulation)) {
      #Passa uma simulação para uma coluna, a abundância já irá em RSA
      matrizLocais[1:length(simLocais$abundance[simLocais$simulation == i]), i] <- simLocais$abundance[simLocais$simulation == i] / sum(simLocais$abundance[simLocais$simulation == i])
    }
    for (x in 1:ncol(matrizLocais)) {
      matrizLocais[,x] <- sort(matrizLocais[,x], decreasing = T)
    }

    #agora passar todos os zeros para NA
    matrizLocais[matrizLocais==0] <- NA

    #retirar a média e desvio da RSA de cada simulação
    RSA <- rowMeans(matrizLocais, na.rm = T)
    SD <- apply(matrizLocais, 1, sd, na.rm=T)
    final <- cbind(RSA, SD)
    final <- as.data.frame(final)
    #agora retirar as linhas que no qual não possua valor de RSA
    final <- final %>% filter(RSA != 0)

    #criando rank, sigma, migração,
    final$RANK <- 1:nrow(final)
    final$SIGMA <- rep(df$sigma[1], nrow(final))
    final$MIG <- rep(df$migration[1], nrow(final))
    print("Método de retirar RSA antes da média aplicado!")
    return(final)
  }

  if (meanRSA == FALSE) {
    for (i in 1: max(simLocais$simulation)) {
      matrizLocais[1:length(simLocais$abundance[simLocais$simulation == i]), i] <- simLocais$abundance[simLocais$simulation == i]
    }
    for (x in 1:ncol(matrizLocais)) {
      matrizLocais[,x] <- sort(matrizLocais[,x], decreasing = T)
    }

    #agora passar todos os zeros para NA
    matrizLocais[matrizLocais==0] <- NA

    #retirar a média e desvio da RSA de cada simulação
    media <- rowMeans(matrizLocais, na.rm = T)
    RSA <- media / sum(media, na.rm = T)
    desvio <- apply(matrizLocais, 1, sd, na.rm=T)
    SD <- desvio/ sum(desvio, na.rm = T)
    final <- cbind(RSA, SD)
    final <- as.data.frame(final)
    #agora retirar as linhas que no qual não possua valor de RSA
    final <- final %>% filter(RSA != 0)

    #criando rank, sigma, migração,
    final$RANK <- 1:nrow(final)
    final$SIGMA <- rep(df$sigma[1], nrow(final))
    final$MIG <- rep(df$migration[1], nrow(final))
    print("Método de retirar a média antes da RSA aplicado!")
    return(final)

  }

}



