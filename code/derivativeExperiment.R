setwd("~/expmDerivative/")

n       <- 100
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)

E        <- matrix(0,n,n)
E[3,1]   <- 1 
frchtDrv <- expm::expmFrechet(A=Q,
                              E=E,
                              expm=FALSE)[[1]]
frchtDrv

frstRdr <- E %*% t(expm::expm(Q)) # correct
norm(frstRdr-frchtDrv)

thrFrstRdr <- expm::expm(Q) %*% E # incorrect
norm(thrFrstRdr-frchtDrv)
