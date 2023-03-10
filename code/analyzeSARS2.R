setwd("~/expmDerivative/")
library(ggplot2)
library(readr)
library(ggthemes)
library(scales)


#
####
####### get predictors
####
#
df <- read_table2("data/airTravel.txt", col_names = FALSE)
airTravs <- unlist(df)
atMat <- matrix(0,44,44)
k <- 1
for(i in 1:43) {
  for(j in (i+1):44) {
    atMat[i,j] <- airTravs[k]
    k <- k+1
  }
}
for(i in 1:43) {
  for(j in (i+1):44) {
    atMat[j,i] <- airTravs[k]
    k <- k+1
  }
}

df <- read_table2("data/continentalDists.txt", col_names = FALSE)
airTravs <- unlist(df)
contMat <- matrix(0,44,44)
k <- 1
for(i in 1:43) {
  for(j in (i+1):44) {
      contMat[i,j] <- airTravs[k]
      k <- k+1
  }
}
for(i in 1:43) {
  for(j in (i+1):44) {
    contMat[j,i] <- airTravs[k]
    k <- k+1
  }
}

df <- read_table2("data/hubeiAsymm.txt", col_names = FALSE)
airTravs <- unlist(df)
hubMat <- matrix(0,44,44)
k <- 1
for(i in 1:43) {
  for(j in (i+1):44) {
    hubMat[i,j] <- airTravs[k]
    k <- k+1
  }
}
for(i in 1:43) {
  for(j in (i+1):44) {
    hubMat[j,i] <- airTravs[k]
    k <- k+1
  }
}

#
###
###### get locations 
###
#
locations <- read_table2("data/locations.txt", col_names = FALSE)
locations[[1]][6] <- "Anhui,CN"
locations[[1]][7] <- "Beijing,CN"
locations[[1]][8] <- "Chongqing,CN"
locations[[1]][9] <- "Fujian,CN"
locations[[1]][10] <- "Guangdong,CN"
locations[[1]][11] <- "Henan,CN"
locations[[1]][12] <- "HK,CN"
locations[[1]][13] <- "Hubei,CN"
locations[[1]][14] <- "Hunan,CN"
locations[[1]][15] <- "Jiangsu,CN"
locations[[1]][16] <- "Jiangxi,CN"
locations[[1]][17] <- "Shandong,CN"
locations[[1]][18] <- "Sichuan,CN"
locations[[1]][19] <- "Yunnan,CN"
locations[[1]][20] <- "Zhejiang,CN"


colnames(atMat) <- unlist(locations)
rownames(atMat) <- unlist(locations)

colnames(contMat) <- unlist(locations)
rownames(contMat) <- unlist(locations)

colnames(hubMat) <- unlist(locations)
rownames(hubMat) <- unlist(locations)

#
###
###### visualize fixed-effects
###
#
atMat <- atMat- (max(atMat)-abs(min(atMat)))/2
atMat <- atMat/max(atMat)
diag(atMat) <- 0
df <- reshape2::melt(atMat)
#df$value <- exp(df$value)
colnames(df) <- c("Var1","Var2","Values")
# df$Values <- df$Values - (max(df$Values)-abs(min(df$Values)))/2
# df$Values <- df$Values / max(df$Values)

gg2 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=Values)) + 
  geom_tile() +
  guides(x =  guide_axis(angle = 90)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Air traffic proximity") +
  scale_fill_gradientn(breaks=c(-1,0,1),
                       colours = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))) +
  theme(plot.title = element_text(vjust = -1),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text())
gg2


contMat <- contMat- (max(contMat)-abs(min(contMat)))/2
contMat <- contMat/max(contMat)
diag(contMat) <- 0
df <- reshape2::melt(contMat)
colnames(df) <- c("Var1","Var2","Values")
gg3 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=`Values`)) + 
  geom_tile() +
  guides(x =  guide_axis(angle = 90)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Intracontinental proximity") +
  scale_fill_gradientn(breaks=c(-1,0,1),
                       colours = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))) +
  theme(plot.title = element_text(vjust = -1),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text())
gg3

df <- reshape2::melt(hubMat)
colnames(df) <- c("Var1","Var2","Values")
gg4 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=`Values`)) + 
  geom_tile() +
  guides(x =  guide_axis(angle = 90)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Hubei asymmetric proximity") +
  scale_fill_gradientn(breaks=c(-1,0,1),
    colours = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))) +
  theme(plot.title = element_text(vjust = -1),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text())
gg4

library(ggpubr)
ggsave(ggpubr::ggarrange(gg2,NULL,gg3,NULL,gg4,NULL,nrow = 1,
                         common.legend=TRUE,
                         legend="right",
                         labels=c("a","","b","","c",""),
                         widths = c(1, -0.02, 1,-0.02,1,-0.02)),
       width = 12,height = 4,
       file="figures/fixedEffects.pdf")
system2(command = "pdfcrop",
        args    = c("~/expmDerivative/figures/fixedEffects.pdf",
                    "~/expmDerivative/figures/fixedEffects.pdf")
)


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
  for(i in 1:43) {
    for(j in (i+1):44) {
        rcMat[i,j] <- unlist(df2[n,k])
        k <- k+1
    }
  }
  for(i in 1:43) {
    for(j in (i+1):44) {
      rcMat[j,i] <- unlist(df2[n,k])
      k <- k+1
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

tC <- reshape2::melt(totalCoeffs[[n-5]])
gg2 <- ggplot(data = tC, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
gg2
  