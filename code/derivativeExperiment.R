setwd("~/expmDerivative/")

#
###
###### correct approximation to frechet derivative
###
#
n       <- 3
Q       <- matrix(rexp(n^2),n,n)#rexp(n^2)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)

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
