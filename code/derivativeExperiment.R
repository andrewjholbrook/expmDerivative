setwd("~/expmDerivative/")

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

# upper bound
ub <- function(t, Q, Jij) {
  n <- dim(Q)[1]
  eig.obj <- eigen(Q)
  inv.spctrm <- 1/eig.obj$values
  inv.spctrm[length(inv.spctrm)] <- 0
  pinv    <- Re(eig.obj$vectors %*% diag(inv.spctrm) %*% solve(eig.obj$vectors))
  mat <- expm::expm(t*Q)%*%Jij%*%pinv%*%(Q+diag(n))
  return(norm(mat,type="F"))
}

nmrclDrvtv <- function(t,Q,E,eps=0.0001) {
  out <- (expm::expm(t*(Q+eps*E/2))-expm::expm(t*(Q-eps*E/2)))/eps
  return(out)
}

n       <- 200
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,3]   <- 1 

ub(1,Q,E)

drv <- nmrclDrvtv(t=1,Q=Q,E=E,eps=0.00001)

frstRdr <- expm::expm(1*Q) %*% E * 1 # correct approx
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
