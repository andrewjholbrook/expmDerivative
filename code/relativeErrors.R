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
n       <- 10
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {

  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx

  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg <- ggplot(df,aes(x=t,y=Error,color=Type)) +
#  geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- n*2
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg2 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  #geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg2

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- n*2
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg3 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  #geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg3


df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- n*2
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg4 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  #geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg4

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- n*2
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg5 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  #geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg5

df      <- data.frame(Error=numeric(),Type=character(),t=numeric())
n       <- n*2
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1

for (t in seq(from=10^(-10),to=10,length.out=100)) {
  
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
  
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F"),
                "Error",
                t))
  df <- rbind(df,
              c(norm(drv,type = "F"),
                "GradNorm",
                t))
  df <- rbind(df,
              c(norm(frstRdr-drv,type = "F")/norm(drv,type = "F"),
                "Relative error",
                t))
  df <- rbind(df,
              c(sqrt(n),
                "sqrt(D)",
                t))
}

colnames(df) <- c("Error","Type","t")
df$Type <- factor(df$Type)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg6 <- ggplot(df,aes(x=t,y=Error,color=Type)) +
  #geom_hline(yintercept = sqrt(n)) +
  geom_line(size=1.5) +
  # scale_y_log10() +
  coord_cartesian(ylim=c(0,18)) +
  ggtitle(paste0("Dimension = ",n)) +
  theme_bw()
gg6

library(ggpubr)
ggsave(ggarrange(gg,gg2,gg3,gg4,gg5,gg6,
                 ncol=2,nrow=3,
                 common.legend = TRUE,
                 legend="bottom"),
       filename = "figures/relativeErrors.pdf",
       device = "pdf",
       width=12,height=12)


