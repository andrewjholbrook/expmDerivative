setwd("~/expmDerivative/")
library(coda)

df <- read_table2("output/short.Ed.TOOBIG4GIT.log", skip = 3)
df <- df[df$state>=1000 & df$state <= 80000,]

# fixed effects
effectiveSize(df[,9:11]) #        1053.3909        1342.5106         498.3175 

# random effects
summary(effectiveSize(df[,12:1903]))
