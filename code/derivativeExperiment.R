setwd("~/expmDerivative/")
library(ggplot2)

# # when do we get convergence?
# n       <- 3
# Q       <- matrix(rexp(n^2),n,n)
# diag(Q) <- 0
# diag(Q) <- - rowSums(Q)
# 
# drv <- expm::expm(10*Q) #stationarity
# drv%*%Q
# 
# # pseudoinverse
# eig.obj <- eigen(Q)
# inv.spctrm <- 1/eig.obj$values
# inv.spctrm[length(inv.spctrm)] <- 0
# pinv    <- eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors)
# norm(pinv%*%Q%*%Q-Q) <= 0.00000001
# 
# drv%*%Q
# 
# drv%*%pinv
# 
# norm(pinv%*%Q,type="F")

# upper bound TODO: show a lower bound
ub <- function(t, Q, Jij) {
  n <- dim(Q)[1]
  eig.obj <- eigen(Q)
  inv.spctrm <- 1/eig.obj$values
  inv.spctrm[length(inv.spctrm)] <- 0
  pinv    <- Re(eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors))
  mat <- expm::expm(t*Q)%*%Jij%*%pinv%*%(Q*t+diag(n))
  return(norm(mat,type="F"))
}

nmrclDrvtv <- function(t,Q,E,eps=0.0001) {
  out <- (expm::expm(t*(Q+eps*E/2))-expm::expm(t*(Q-eps*E/2)))/eps
  return(out)
}

blockwise <- function(t,Q,E) {
  B1  <- cbind(Q,E)
  B2  <- cbind(Q*0,Q)
  B   <- rbind(B1,B2) 
  etB <- expm::expm(t*B)
  return(etB[1:(dim(etB)[1]/2),(dim(etB)[1]/2+1):dim(etB)[1]])
}

# 
######## figure
#
n       <- 50
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
eigs    <- eigen(Q)$values

for(n in 100*1:10) {
  Q       <- matrix(rexp(n^2),n,n)
  diag(Q) <- 0
  diag(Q) <- - rowSums(Q)
  eigs    <- c(eigs,eigen(Q)$values)
}

df <- data.frame(Real=Re(eigs),
                 Imaginary=Im(eigs),
                 Dimension=c(rep(50,50),
                             rep(100,100),
                             rep(200,200),
                             rep(300,300),
                             rep(400,400),
                             rep(500,500),
                             rep(600,600),
                             rep(700,700),
                             rep(800,800),
                             rep(900,900),
                             rep(1000,1000)))
df$Dimension <- factor(df$Dimension)

gg <- ggplot(df,aes(x=Real,y=Imaginary,color=Dimension)) +
  geom_point(size=0.1) +
  ggtitle("Eigenvalues of random CTMC generators") +
  theme_bw()
gg

ggsave("eigenvaluesPlot.png",device = "png",width=5,height=3, dpi=320)

##########################

n       <- 100
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,2]   <- 1 
t <- 1000000000

#max(abs(expm::expm(t*Q)%*%Q))

ub(t,Q,E)

#drv <- nmrclDrvtv(t=t,Q=Q,E=E,eps=0.00001)
drv <- blockwise(t=t,Q,E)

frstRdr <- expm::expm(t*Q) %*% E *t # correct approx
norm(frstRdr-drv,type = "F")


#
###
####### concentration of measure?
###
#
n       <- 10000
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
sum(abs(Q[2:n,1]))

head(sort(Mod(eigen(Q)$values)))

n       <- 1000
Q       <- matrix(exp(rnorm(n=n^2,mean=0,sd=10)),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
head(sort(Mod(eigen(Q)$values)))

#
###
###### correct approximation to frechet derivative
###
#
n       <- 3
Q       <- matrix(rexp(n^2),n,n)#rexp(n^2)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)


############

E        <- matrix(0,n,n)
E[3,2]   <- 1 
frchtDrv <- expm::expmFrechet(A=t(Q),
                              E=E,
                              expm=FALSE)[[1]]

frstRdr <- E %*% t(expm::expm(Q)) # correct approx
frstRdr
norm(frstRdr-frchtDrv,type = "F")
max(abs(frstRdr-frchtDrv))

thrFrstRdr <-Q%*% expm::expm(1*Q) %*% E # incorrect(?) approx
thrFrstRdr[,2]
sd(thrFrstRdr[,1])
norm(thrFrstRdr-frchtDrv,type = "F")
max(abs(thrFrstRdr-frchtDrv))
#
####
####### compare frechet derivative to numerical gradient
####
#

nmrclDrvtv <- function(Q,E,eps=0.0001) {
  out <- (expm::expm(Q+eps*E/2)-expm::expm(Q-eps*E/2))/eps
  return(out)
}

n       <- 100
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)

E        <- matrix(0,n,n)
E[3,1]   <- 1 

nmDrv <- nmrclDrvtv(Q,E,eps=0.0001)

# correct: pertains to perturbing single entry
thrFrchtDrv <- expm::expmFrechet(A=Q,
                              E=E,
                              expm=FALSE)[[1]]
norm(thrFrchtDrv-nmDrv,type="F")
max(abs(thrFrchtDrv-nmDrv))


# incorrect?
frstRdr <- E %*% t(expm::expm(Q)) 
norm(frstRdr-thrFrchtDrv,type="F")
max(abs(frstRdr-thrFrchtDrv))


# correct
thrFrstRdr <- expm::expm(Q) %*% E 
norm(thrFrstRdr-thrFrchtDrv,type="F")
max(abs(thrFrstRdr-thrFrchtDrv))



#
###
###### structured directional derivative
###
#
n       <- 3
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[3,1]   <- 1 
E[3,3]   <- -1
frchtDrv <- expm::expmFrechet(A=Q,
                              E=E,
                              expm=FALSE)[[1]]

frstRdr <- expm::expm(Q) %*% E # correct
norm(frstRdr-frchtDrv,type="F") # should be 0???
