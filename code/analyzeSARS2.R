setwd("~/expmDerivative/")
library(ggplot2)
library(readr)

#
####
####### get predictors
####
#
df <- read_table2("data/airTravel.txt", col_names = FALSE)
airTravs <- unlist(df)
atMat <- matrix(0,44,44)
k <- 1
for(i in 1:44) {
  for(j in 1:44) {
    if(i!=j){
      atMat[i,j] <- airTravs[k]
      k <- k+1
    }
  }
}

df <- read_table2("data/continentalDists.txt", col_names = FALSE)
airTravs <- unlist(df)
contMat <- matrix(0,44,44)
k <- 1
for(i in 1:44) {
  for(j in 1:44) {
    if(i!=j){
      contMat[i,j] <- airTravs[k]
      k <- k+1
    }
  }
}

df <- read_table2("data/hubeiAsymm.txt", col_names = FALSE)
airTravs <- unlist(df)
hubMat <- matrix(0,44,44)
k <- 1
for(i in 1:44) {
  for(j in 1:44) {
    if(i!=j){
      hubMat[i,j] <- airTravs[k]
      k <- k+1
    }
  }
}

#
###
###### get locations (incorrect)
###
#
locations <- read_table2("data/locations.txt", col_names = FALSE)

colnames(atMat) <- unlist(locations)
rownames(atMat) <- unlist(locations)

colnames(contMat) <- unlist(locations)
rownames(contMat) <- unlist(locations)

colnames(hubMat) <- unlist(locations)
rownames(hubMat) <- unlist(locations)

#
###
###### visualize posts
###
#
df <- read_table2("output/covid10Mar_285plusNYC.Ed.log", skip = 3)

df2 <- df[,9:11]
df2 <- reshape2::melt(df2)
df2$Predictor <- "cat"
df2$Predictor[which(df2$variable=="glmCoefficients1")] <- "IntraContinental"
df2$Predictor[which(df2$variable=="glmCoefficients2")] <- "AirTraffic"
df2$Predictor[which(df2$variable=="glmCoefficients3")] <- "HubeiAsymmetry"

gg <- ggplot(df2,aes(x=value,color=Predictor)) +
  geom_density(adjust=2) +
  xlab("Coefficient value") +
  ylab("Density") +
  ggtitle("Fixed-effects for infinitesimal rates between locations") +
  theme_bw()
gg

#
####
####### get random coeffs matrices
####
#
df2 <- df[,12:1903]
nSamps <- dim(df2)[1]
randCoeffs <- list()


for(n in 1:nSamps) {
  k <- 1
  rcMat <- matrix(0,44,44)
  for(i in 1:44) {
    for(j in 1:44) {
      if(i!=j){
        rcMat[i,j] <- unlist(df2[n,k])
        k <- k+1
      }
    }
  }
  randCoeffs[[n]] <- rcMat
}

totalCoeffs <- list()
for (n in 1:nSamps) {
  totCoeffMat <- randCoeffs[[n]] +
    atMat * unlist(df[n,10]) +
    hubMat * unlist(df[n,11]) +
    contMat * unlist(df[n,9])
  
  totalCoeffs[[n]] <- exp(totCoeffMat)
  diag(totalCoeffs[[n]]) <- 0
}
  
  
  