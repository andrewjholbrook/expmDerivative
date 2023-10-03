setwd("~/expmDerivative/")

##### simulation design (3x3x3)
# 1) 3 significance proportions: 5 percent, 50 percent, 95 percent
# 2) 3 gradients: exact, affine-corrected approx, naive approx
# 3) 3 state space sizes: 15, 30, 45

library(expm)

################################################################################
# code for simulating data

simCTMC <- function(Q, initialStates, time) {
  S <- dim(Q)[1]
  N <- length(initialStates)
  if(any(initialStates>S)) stop("Initial states cannot be larger than state space.")
  
  transMat <- expm(time*Q)
  
  output <- rep(0,N)
  for(n in 1:N) {
    ntlStt <- rep(0,S)
    ntlStt[initialStates[n]] <- 1
    outputProbs <- t(transMat) %*% ntlStt
    output[n] <- sample(x=1:S,size=1,prob=outputProbs)
  }
  return(output)
}

# # test
# S <- 15
# Q       <- matrix(rexp(S^2),S,S)
# diag(Q) <- 0
# diag(Q) <- - rowSums(Q)
# initialStates <- sample(x=1:S,size=100,replace=TRUE)
# simCTMC(Q,initialStates,1)

lg_lklhd <- function(Q,input,output,time) {
  if(! length(input)==length(output)) stop("Input and output must have same length.")
  S <- dim(Q)[1]
  if(any(input>S)) stop("Input states cannot be larger than state space.")
  if(any(output>S)) stop("Output states cannot be larger than state space.")
  N <- length(input)
  transMat <- expm::expm(time*Q)
  
  lklhd <- 0
  for(n in 1:N) {
    ntlStt <- rep(0,S)
    tptStt <- rep(0,S)
    ntlStt[input[n]] <- 1
    tptStt[output[n]] <- 1
    lklhd <- lklhd + log(t(tptStt) %*% t(transMat) %*% ntlStt)
  }
  return(as.numeric(lklhd))
}

# # test
# S <- 15
# Q       <- matrix(rexp(S^2),S,S)
# diag(Q) <- 0
# diag(Q) <- - rowSums(Q)
# initialStates <- sample(x=1:S,size=100,replace=TRUE)
# output <- simCTMC(Q,initialStates,1)
# 
# ll <- lg_lklhd(Q,initialStates,output,1)
# ll

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

frstRdr <- function(t,Q,E){
  return(expm::expm(t*Q) %*% E * t)
}

frstRdrCrrctd <- function(t,Q,E) {
  out <- frstRdr(t,Q,E) - approx(t,Q,E)
  return(out)
}

lg_prr <- function(Q) {
  return(0) # TODO: add bayesian bridge log prior
}

lg_prr_grd <- function(Q) {
  return(0) # TODO: add bayesian bridge log prior gradient
}

lg_lklhd_grd <- function(t,Q,input,output,method) {
  if(! length(input)==length(output)) stop("Input and output must have same length.")
  S <- dim(Q)[1]
  if(any(input>S)) stop("Input states cannot be larger than state space.")
  if(any(output>S)) stop("Output states cannot be larger than state space.")
  N <- length(input)
  if(! method %in% c("exact","approx","corrected")) stop("Use prescribed gradient method!")
  gradient <- matrix(0,S,S)
  
  for(i in 1:S) {
    Eii <- matrix(0,S,S)
    Eii[i,i] <- 1
    for(j in 1:S) {
      if(i != j) {
        patty <- matrix(0,S,S)
        E <- matrix(0,S,S)
        E[i,j] <- 1
        if (method=="exact") {
          patty <- patty + blockwise(t,Q,E) * Q[i,j] + blockwise(t,Q,Eii) * Q[i,i]
        } else if (method=="approx") {
          patty <- patty + frstRdr(t,Q,E) * Q[i,j] + frstRdr(t,Q,Eii) * Q[i,i]
        } else {
          patty <- patty + frstRdrCrrctd(t,Q,E) * Q[i,j] + frstRdrCrrctd(t,Q,Eii) * Q[i,i]
        }
        for(n in 1:N) {
          ntlStt <- rep(0,S)
          tptStt <- rep(0,S)
          ntlStt[input[n]] <- 1
          tptStt[output[n]] <- 1
          gradient[i,j] <- gradient[i,j] + t(tptStt) %*% t(patty) %*% ntlStt
        }
      }
    }
  }
  
  # don't forget to divide by log-likelihood
  gradient <- gradient / lg_lklhd(Q,input,output,t) 
  diag(gradient) <- 0
  
  return(gradient)
 }

# test
# S <- 5
# Q       <- matrix(rexp(S^2),S,S)
# diag(Q) <- 0
# diag(Q) <- - rowSums(Q)
# initialStates <- sample(x=1:S,size=100,replace=TRUE)
# output <- simCTMC(Q,initialStates,1)
# 
# ll_grd <- lg_lklhd_grd(1,Q,initialStates,output,method="exact")
# ll_grd
# 
# ll_grd_prx <- lg_lklhd_grd(1,Q,initialStates,output,method="approx")
# ll_grd_prx
# 
# ll_grd_crctd <- lg_lklhd_grd(1,Q,initialStates,output,method="corrected")
# ll_grd_crctd
# 
# norm(ll_grd - ll_grd_prx)
# 
# norm(ll_grd - ll_grd_crctd)
# 
# S <- 45
# Q       <- matrix(rexp(S^2),S,S)
# diag(Q) <- 0
# diag(Q) <- - rowSums(Q)
# initialStates <- sample(x=1:S,size=100,replace=TRUE)
# output <- simCTMC(Q,initialStates,1)
# 
# ll_grd <- lg_lklhd_grd(1,Q,initialStates,output,method="exact")
# #ll_grd
# 
# ll_grd_prx <- lg_lklhd_grd(1,Q,initialStates,output,method="approx")
# #ll_grd_prx
# 
# ll_grd_crctd <- lg_lklhd_grd(1,Q,initialStates,output,method="corrected")
# #ll_grd_crctd
# 
# norm(ll_grd - ll_grd_prx)
# 
# norm(ll_grd - ll_grd_crctd)

##########################################################
# HMC
target <- function(Q,input,output,time) {
  out <- lg_lklhd(Q,input,output,time) + lg_prr(Q)
  return(out)
}

grad <- function(t,Q,input,output,method) {
  out <- lg_lklhd_grd(t,Q,input,output,method) + lg_prr_grd(Q)
  return(out)
}

delta <- function(n) {
  #return( min(0.01,n^(-0.5)) )
  return(n^(-0.5))
}

hmc <- function(t,S,input,output,method,stepSize=0.01,L=20,maxIts=1000,
                targetAccept=0.8) {
  
  Qs <- list()
  Qs[[1]] <- matrix(rexp(S^2),S,S)
  diag(Qs[[1]]) <- 0
  diag(Qs[[1]]) <- - rowSums(Qs[[1]])

  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  
  currentU  <- - target(Qs[[1]],input,output,t)
  log_probs <- rep(0,maxIts)
  log_probs[1] <- currentU
  
  for (i in 2:maxIts) {
    proposalState    <- Qs[[i-1]]
    momentum         <- matrix(rnorm(S^2),S,S)
    diag(momentum)   <- 0
    currentK         <- sum(momentum^2)/2
    
    # leapfrog steps
    momentum <- momentum + 0.5 * stepSize * grad(t,proposalState,input,output,method)
    for (l in 1:L) {
      proposalState <- proposalState + stepSize * momentum
      diag(proposalState) <- 0
      diag(proposalState) <- - rowSums(proposalState) 
      if (l!=L) momentum <- momentum + stepSize * grad(t,proposalState,input,output,method)
    }
    momentum <- momentum + 0.5 * stepSize * grad(t,proposalState,input,output,method)
    
    # quantities for accept/reject
    proposedU = - target(proposalState,input,output,t)
    proposedK = sum(momentum^2)/2
    u <- runif(1)
    
    if (log(u) < currentU + currentK - proposedU - proposedK) {
      Qs[[i]]     <- proposalState
      currentU    <- proposedU
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    } else {
      Qs[[i]] <- Qs[[i-1]]
    }
    
    SampCount <- SampCount + 1
    log_probs[i] <- currentU
    
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        stepSize <- stepSize * (1 + delta(i-1))
      } else {
        stepSize <- stepSize * (1 - delta(i-1))
      }
      
      SampCount <- 0
      Acceptances <- 0
    }
    
    
    if (i %% 100 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n") 
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(list(Qs,-log_probs))
}

# test
S <- 5
maxIts <- 1000
initialStates <- sample(x=1:S,size=100,replace=TRUE)
Q       <- matrix(rexp(S^2),S,S)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
output <- simCTMC(Q,initialStates,1)
hmcOut <- hmc(t=1,S=S,input=initialStates,output=output,stepSize=0.005,
             method="exact",maxIts = maxIts, targetAccept = 0.65)
plot(hmcOut[[2]],type="l")
chain <- hmcOut[[1]]

variable <- rep(0,maxIts)
for(i in 1:maxIts) {
  variable[i] <- log(chain[[i]][2,1])
}
plot(variable,type="l")
abline(h=log(Q[2,1]),col="red")

