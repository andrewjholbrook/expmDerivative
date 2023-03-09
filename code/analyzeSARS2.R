setwd("~/expmDerivative/")
library(ggplot2)
library(readr)

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

