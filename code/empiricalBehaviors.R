setwd("~/expmDerivative/")
library(ggplot2)

blockwise <- function(t,Q,E) {
  B1  <- cbind(Q,E)
  B2  <- cbind(Q*0,Q)
  B   <- rbind(B1,B2) 
  etB <- expm::expm(t*B)
  return(etB[1:(dim(etB)[1]/2),(dim(etB)[1]/2+1):dim(etB)[1]])
}

cnvrgncBnd <- function(t, Q, Jij) {
  n <- dim(Q)[1]
  eig.obj <- eigen(Q)
  inv.spctrm <- 1/eig.obj$values
  inv.spctrm[length(inv.spctrm)] <- 0
  pinv    <- Re(eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors))
  mat <- expm::expm(t*Q)%*%Jij%*%pinv%*%(Q*t+diag(n))
  return(norm(mat,type="F"))
}

nvBnd1 <- function(t, Q, Jij) {
  nrmXptQ <- norm(expm::expm(t*Q),type="F")
  nrmQ    <- norm(Q,type="F")
  out     <- nrmXptQ / (2*nrmQ) * (exp(2*t*nrmQ) - 2*t*nrmQ - 1)
  return(out)
}

nvBnd2 <- function(t, Q, Jij) {
  nrmQ    <- norm(Q,type="F")
  out     <- sqrt(2)/nrmQ * (exp(t*nrmQ)*(t*nrmQ - 1) + 1)
  return(out)
}

################################################################################


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 100
Q       <- matrix(1,n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 
vals <- diag(1/eigen(Q)$values)
vals[1,1] <- 0
vecs <- eigen(Q)$vectors
#max(abs(vecs%*%vals%*%solve(vecs)-Q))
Qplus <- vecs%*%vals%*%solve(vecs)
norm(Qplus)

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  geom_hline(yintercept = 10*norm(Qplus,type="2")) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 5") +
  theme_bw()
gg

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 10
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[1,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg2 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
 scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 10") +
  theme_bw()
gg2


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 15
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg3 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 15") +
  theme_bw()
gg3

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 20
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg4 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 20") +
  theme_bw()
gg4

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 25
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg5 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 25") +
  theme_bw()
gg5


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 30
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg6 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 30") +
  theme_bw()
gg6

library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,gg4,gg5,gg6,
                 ncol=2,nrow=3,
                 common.legend = TRUE,
                 legend="bottom"),
       filename = "frobeniusScaling.pdf",
       device = "pdf",
       width=12,height=12)


################################################################################
# replace with i norm

blockwise <- function(t,Q,E) {
  B1  <- cbind(Q,E)
  B2  <- cbind(Q*0,Q)
  B   <- rbind(B1,B2) 
  etB <- expm::expm(t*B)
  return(etB[1:(dim(etB)[1]/2),(dim(etB)[1]/2+1):dim(etB)[1]])
}

cnvrgncBnd <- function(t, Q, Jij) {
  n <- dim(Q)[1]
  eig.obj <- eigen(Q)
  inv.spctrm <- 1/eig.obj$values
  inv.spctrm[length(inv.spctrm)] <- 0
  pinv    <- Re(eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors))
  mat <- expm::expm(t*Q)%*%Jij%*%pinv%*%(Q*t+diag(n))
  return(norm(mat,type="i"))
}

nvBnd1 <- function(t, Q, Jij) {
  nrmXptQ <- norm(expm::expm(t*Q),type="i")
  nrmQ    <- norm(Q,type="F")
  out     <- nrmXptQ / (2*nrmQ) * (exp(2*t*nrmQ) - 2*t*nrmQ - 1)
  return(out)
}

nvBnd2 <- function(t, Q, Jij) {
  nrmQ    <- norm(Q,type="i")
  out     <- sqrt(2)/nrmQ * (exp(t*nrmQ)*(t*nrmQ - 1) + 1)
  return(out)
}


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 5
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 5") +
  theme_bw()
gg

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 10
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg2 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 10") +
  theme_bw()
gg2


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 15
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg3 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 15") +
  theme_bw()
gg3

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 20
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg4 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 20") +
  theme_bw()
gg4

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 25
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg5 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 25") +
  theme_bw()
gg5


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- 30
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1 

for (t in seq(from=10^(-10),to=1,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Truth",
                t))
  df <- rbind(df,
              c(cnvrgncBnd(t,Q,E),
                "CnvgAssumption",
                t))
  df <- rbind(df,
              c(nvBnd1(t,Q,E),
                "Bound1",
                t))
  df <- rbind(df,
              c(nvBnd2(t,Q,E),
                "Bound2",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg6 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  geom_line(size=1.5) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-4,1e+4)) +
  ggtitle("Dimension = 30") +
  theme_bw()
gg6

library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,gg4,gg5,gg6,
                 ncol=2,nrow=3,
                 common.legend = TRUE,
                 legend="bottom"),
       filename = "infinityScaling.pdf",
       device = "pdf",
       width=12,height=12)
