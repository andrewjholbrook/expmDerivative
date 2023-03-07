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

df      <- data.frame(Error=numeric(),Dimension=character(),t=numeric())


for(n in 3:50) {
for (t in seq(from=10^(-10),to=10,length.out=100)) {

  Q       <- matrix(rexp(n^2),n,n)
  diag(Q) <- 0
  diag(Q) <- - rowSums(Q)
  E       <- matrix(0,n,n)
  E[sample(1:n,size=1),sample(1:n,size=1)]   <- 1
  #d <- sample(1:n,size=1)
  #E[d,d] <- 1
  drv <- blockwise(t=t,Q,E)
  frstRdr <- expm::expm(t*Q) %*% E *t # correct approx

  df <- rbind(df,
              c(acos(sum(diag(t(frstRdr)%*%drv))/
                  (norm(frstRdr,type="F")*norm(drv,type="F")))*(360/(2*pi)),
                n,
                t))

}
}

colnames(df) <- c("Error","Dimension","t")
df$Dimension <- factor(df$Dimension)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg <- ggplot(df,aes(x=t,y=Error,color=Dimension)) +
  geom_line(size=1.5) +
  geom_hline(yintercept = 90) +
  # scale_y_log10() +
  ylab("Angle (deg.)") +
  coord_cartesian(ylim=c(0,100)) +
  ggtitle("Directional derivative vs 1st order approx.") +
  theme_bw()
gg

# df      <- data.frame(Error=numeric(),Dimension=character(),t=numeric())
# 
# 
# for(n in 3:50) {
#   for (t in seq(from=10^(-10),to=10,length.out=100)) {
#     
#     Q       <- matrix(rexp(n^2),n,n)
#     diag(Q) <- 0
#     diag(Q) <- - rowSums(Q)
#     E       <- matrix(0,n,n)
#     d <- sample(1:n,size=1)
#     dprime <- sample((1:n)[-d],size=1)
#     E[d,dprime]   <- 1
#     #E[2,2] <- 1
#     drv <- blockwise(t=t,Q,E)
#     frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
#     
#     df <- rbind(df,
#                 c(acos(sum(diag(t(frstRdr)%*%drv))/
#                          (norm(frstRdr,type="F")*norm(drv,type="F")))*(360/(2*pi)),
#                   n,
#                   t))
#     
#   }
# }
# 
# colnames(df) <- c("Error","Dimension","t")
# df$Dimension <- factor(df$Dimension)
# df$Error <- as.numeric(df$Error)
# df$t     <- as.numeric(df$t)
# 
# gg2 <- ggplot(df,aes(x=t,y=Error,color=Dimension)) +
#   geom_line(size=1.5) +
#   geom_hline(yintercept = 90) +
#   # scale_y_log10() +
#   ylab(" ") +
#   coord_cartesian(ylim=c(0,100)) +
#   ggtitle("Off-diagonal directional derivative vs 1st order approx.") +
#   theme_bw()
# gg2

df      <- data.frame(Error=numeric(),Dimension=character(),t=numeric())


for(n in 3:50) {
  for (t in seq(from=10^(-10),to=10,length.out=100)) {

    Q       <- matrix(rexp(n^2),n,n)
    diag(Q) <- 0
    diag(Q) <- - rowSums(Q)
    E       <- matrix(0,n,n)
    d <- sample(1:n,size=1)
    dprime <- sample((1:n)[-d],size=1)
    E[d,dprime]   <- 1
    E[d,d]  <- -1
    #E[2,2] <- 1
    drv <- blockwise(t=t,Q,E)
    frstRdr <- expm::expm(t*Q) %*% E *t # correct approx

    df <- rbind(df,
                c(acos(sum(diag(t(frstRdr)%*%drv))/
                         (norm(frstRdr,type="F")*norm(drv,type="F")))*(360/(2*pi)),
                  n,
                  t))

  }
}

colnames(df) <- c("Error","Dimension","t")
df$Dimension <- factor(df$Dimension)
df$Error <- as.numeric(df$Error)
df$t     <- as.numeric(df$t)

gg3 <- ggplot(df,aes(x=t,y=Error,color=Dimension)) +
  geom_line(size=1.5) +
  geom_hline(yintercept = 90) +
  # scale_y_log10() +
  ylab(" ") +
  coord_cartesian(ylim=c(0,100)) +
  ggtitle("\"Normalized\" directional derivative vs 1st order approx.") +
  theme_bw()
gg3

library(ggpubr)
ggsave(ggarrange(gg,gg3,
                 ncol=2,nrow=1,
                 common.legend = TRUE,
                 legend="right"),
       filename = "figures/angles.pdf",
       device = "pdf",
       width=12,height=4)


