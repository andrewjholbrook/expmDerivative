setwd("~/expmDerivative/")

n       <- 2000
Q       <- matrix(rexp(n^2),n,n)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
E        <- matrix(0,n,n)
E[2,1]   <- 1
vals <- diag(1/eigen(Q)$values)
vals[1,1] <- 0
vecs <- eigen(Q)$vectors
#max(abs(vecs%*%vals%*%solve(vecs)-Q))
Qplus <- vecs%*%vals%*%solve(vecs)
max(Mod(Q%*%Qplus-diag(n)))
