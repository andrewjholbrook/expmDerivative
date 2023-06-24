setwd("~/expmDerivative/")
library(ggplot2)
library(ggthemes)

blockwise <- function(t,Q,E) {
  B1  <- cbind(Q,E)
  B2  <- cbind(Q*0,Q)
  B   <- rbind(B1,B2) 
  etB <- expm::expm(t*B)
  return(etB[1:(dim(etB)[1]/2),(dim(etB)[1]/2+1):dim(etB)[1]])
}

approx <- function(t, Q, Jij) {
  n <- dim(Q)[1]
  eig.obj <- eigen(Q)
  inv.spctrm <- 1/eig.obj$values
  inv.spctrm[length(inv.spctrm)] <- 0
  pinv    <- Re(eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors))
  mat <- t*(diag(n) - pinv%*%Q) %*% Jij %*% Q %*% pinv
  return(mat)
}

################################################################################
set.seed(1)
df      <- data.frame(Dimension=numeric(),
                      Error=numeric(),
                      Type=character(),
                      t=numeric())
for(k in 1:20) {
for(n in c(5,20,40,80,160)) {
  Q       <- matrix(rexp(n^2),n,n)
  #Q       <- 0.5* (Q+t(Q))
  diag(Q) <- 0
  diag(Q) <- - rowSums(Q)
  E        <- matrix(0,n,n)
  E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1
  
  for (t in seq(from=10^(-10),to=1,length.out=10)) {
    
    drv <- blockwise(t=t,Q,E)
    frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
    
    df <- rbind(df,
                c(n, norm(frstRdr-drv,type = "F"),
                  "Uncorrected",
                  t))
    df <- rbind(df,
                c(n,norm(frstRdr-approx(t,Q,E)-drv,type = "F"),
                  "Corrected",
                  t))
  }
}
}

colnames(df) <- c("Dimension","Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)
df$Dimension <- factor(df$Dimension,levels=c("5","20","40","80","160"))
#levels(df$Dimension) <- c("5","20","40","80","160")
tail(df)

gg <- ggplot(df,aes(x=t,y=Error,color=Dimension,linetype=Type)) +
  geom_smooth(size=1.5,se=FALSE) + ylab("Approximation error") +
  scale_color_manual(values = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))[c(1,7,13,19,25)])+
  theme_bw() +
  ggtitle("Asymmetric and iid") +
  theme(legend.text = element_text(colour="black"))+
  coord_cartesian(ylim=c(0,0.35))

gg




df      <- data.frame(Dimension=numeric(),
                      Error=numeric(),
                      Type=character(),
                      t=numeric())
for(k in 1:20) {
  for(n in c(5,20,40,80,160)) {
    Q       <- matrix(rexp(n^2),n,n)
    Q       <- Q + matrix(rexp(n),n,n,byrow = TRUE)
    Q       <- Q + matrix(rexp(n),n,n,byrow = FALSE)
    #Q       <- 0.5* (Q+t(Q))
    diag(Q) <- 0
    diag(Q) <- - rowSums(Q)
    E        <- matrix(0,n,n)
    E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1
    
    for (t in seq(from=10^(-10),to=1,length.out=10)) {
      
      drv <- blockwise(t=t,Q,E)
      frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
      
      df <- rbind(df,
                  c(n, norm(frstRdr-drv,type = "F"),
                    "Uncorrected",
                    t))
      df <- rbind(df,
                  c(n,norm(frstRdr-approx(t,Q,E)-drv,type = "F"),
                    "Corrected",
                    t))
    }
  }
}

colnames(df) <- c("Dimension","Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)
df$Dimension <- factor(df$Dimension,levels=c("5","20","40","80","160"))
#levels(df$Dimension) <- c("5","20","40","80","160")
tail(df)

gg2 <- ggplot(df,aes(x=t,y=Error,color=Dimension,linetype=Type)) +
  geom_smooth(size=1.5,se=FALSE) + ylab("") +
  scale_color_manual(values = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))[c(1,7,13,19,25)])+
  theme_bw() +
  ggtitle("Asymmetric and independent") +
  theme(legend.text = element_text(colour="black"))+
  coord_cartesian(ylim=c(0,0.35))
gg2


library(ggpubr)




df      <- data.frame(Dimension=numeric(),
                      Error=numeric(),
                      Type=character(),
                      t=numeric())
for(k in 1:20) {
  for(n in c(5,20,40,80,160)) {
    Q       <- matrix(abs(rt(n=n^2,df=1)),n,n)
    #Q       <- 0.5* (Q+t(Q))
    diag(Q) <- 0
    diag(Q) <- - rowSums(Q)
    E        <- matrix(0,n,n)
    E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1
    
    for (t in seq(from=10^(-10),to=1,length.out=10)) {
      
      drv <- blockwise(t=t,Q,E)
      frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
      
      df <- rbind(df,
                  c(n, norm(frstRdr-drv,type = "F"),
                    "Uncorrected",
                    t))
      df <- rbind(df,
                  c(n,norm(frstRdr-approx(t,Q,E)-drv,type = "F"),
                    "Corrected",
                    t))
    }
  }
}

colnames(df) <- c("Dimension","Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)
df$Dimension <- factor(df$Dimension,levels=c("5","20","40","80","160"))
#levels(df$Dimension) <- c("5","20","40","80","160")
tail(df)

gg3 <- ggplot(df,aes(x=t,y=Error,color=Dimension,linetype=Type)) +
  geom_smooth(size=1.5,se=FALSE) + ylab("") +
  scale_color_manual(values = tableau_div_gradient_pal("Temperature Diverging")(seq(0, 1, length = 25))[c(1,7,13,19,25)])+
  theme_bw() +
  ggtitle("Asymmetric, iid and heavy-tailed") +
  theme(legend.text = element_text(colour="black"))+
  coord_cartesian(ylim=c(0,0.35))
gg3

library(ggpubr)
ggsave(ggpubr::ggarrange(gg,NULL,gg2,NULL,gg3,nrow = 1,
                         common.legend=TRUE,
                         legend="bottom",
                         labels=c("a","","b","","c"),
                         widths = c(1, -0.02, 1,-0.02,1)),
       width = 12,height = 4,
       file="figures/appendixFig3.pdf")
system2(command = "pdfcrop",
        args    = c("~/expmDerivative/figures/appendixFig3.pdf",
                    "~/expmDerivative/figures/appendixFig3.pdf")
)
