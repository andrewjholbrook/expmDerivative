setwd("~/expmDerivative/")

library(readr)
library(ggplot2)

df <- read_table(paste0("data/highDimSim_",10,"_",2,".txt"),
                       col_names = FALSE)
P <- dim(df)[2] - 1
df <- df[,1:P]
df <- df[,unlist(df[1,]) != 0]
Truth <- unlist(df[1,])
df <- df[complete.cases(df),]
N <- dim(df)[1]
df2 <- data.frame(Truth=Truth,
                  Mean=unlist(colMeans(df[4001:5000,])),
                  Dimension=10,
                  Exponent=2)

  
  
for (Dimension in c(10,20,30) ) {
  for(alpha in c(2,1,"half","fourth")) {
    df <- read_table(paste0("data/highDimSim_",Dimension,"_",alpha,".txt"),
                     col_names = FALSE)
    P <- dim(df)[2] - 1
    df <- df[,1:P]
    df <- df[,unlist(df[1,]) != 0]
    Truth <- unlist(df[1,])
    df <- df[complete.cases(df),]
    N <- dim(df)[1]
    
    Exponent <- alpha
    if(Exponent=="half") {
      Exponent <- 0.5
    } else if (Exponent=="fourth") {
      Exponent <- 0.25
    }
    
    df_temp <- data.frame(Truth=Truth,
                      Mean=unlist(colMeans(df[4001:5000,])),
                      Dimension=Dimension,
                      Exponent=Exponent)
    
    df2 <- rbind(df2,df_temp)
  }
}

df <- df2[!duplicated(df2),]
remove(df2)

df$Dimension <- factor(df$Dimension)
df$Exponent  <- factor(df$Exponent)

set.seed(1)
df <- df[sample(1:dim(df)[1]),]

gg <- ggplot(df[df$Exponent==0.25,],aes(x=Truth,y=Mean,color=Dimension)) +
  geom_abline(aes(slope=1,intercept=0),color="black") +
  geom_jitter(size=3) +
  scale_color_manual(values = ggthemes::stata_pal("s2color")(15)[c(7,3,4)]) +
  coord_cartesian(ylim=c(-1.1,1.1),xlim=c(-1.1,1.1)) +
  ylab("") +
  ggtitle("Bayesian bridge exponent = 0.25") +
  theme_bw()
gg

gg2 <- ggplot(df[df$Exponent==0.5,],aes(x=Truth,y=Mean,color=Dimension)) +
  geom_abline(aes(slope=1,intercept=0),color="black") +
  geom_jitter(size=3) +
  scale_color_manual(values = ggthemes::stata_pal("s2color")(15)[c(7,3,4)]) +
  coord_cartesian(ylim=c(-1.1,1.1),xlim=c(-1.1,1.1)) +
  ylab("") +
  ggtitle("Bayesian bridge exponent = 0.5") +
  theme_bw()
gg2

gg3 <- ggplot(df[df$Exponent==1,],aes(x=Truth,y=Mean,color=Dimension)) +
  geom_abline(aes(slope=1,intercept=0),color="black") +
  geom_jitter(size=3) +
  scale_color_manual(values = ggthemes::stata_pal("s2color")(15)[c(7,3,4)]) +
  coord_cartesian(ylim=c(-1.1,1.1),xlim=c(-1.1,1.1)) +
  ylab("Posterior mean") +
  ggtitle("Bayesian bridge exponent = 1") +
  theme_bw()
gg3

library(ggpubr)

ggsave(ggpubr::ggarrange(gg3,NULL,gg2,NULL,gg,nrow = 1,
                         common.legend=TRUE,
                         legend="bottom",
                         labels=c("a","","b","","c"),
                         widths = c(1, -0.02, 1,-0.02,1)),
       width = 12,height = 4,
       file="figures/truthMeanFig.pdf")
system2(command = "pdfcrop",
        args    = c("~/expmDerivative/figures/truthMeanFig.pdf",
                    "~/expmDerivative/figures/truthMeanFig.pdf")
)






