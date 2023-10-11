setwd("~/expmDerivative/")

library(readr)
library(coda)

df <- read_table("data/highDimSim_20_half.txt", col_names = FALSE)
N <- dim(df)[1]
P <- dim(df)[2] - 1
df <- df[,1:P]
df <- df[,as.vector(df[1,]!=0)]
plot(as.vector(df[1,]),as.vector(colMeans(df)))
abline(a=0,b=1,col="red")
cor(as.numeric(df[1,]),as.numeric(colMeans(df)))


plot(1:N,as.numeric(unlist(df[,1])))
