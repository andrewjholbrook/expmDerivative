setwd("~/expmDerivative/")

##### simulation design (3x3x3)
# 1) 3 significance proportions: 5 percent, 50 percent, 95 percent
# 2) 3 gradients: exact, affine-corrected approx, naive approx
# 3) 3 state space sizes: 15, 30, 45

library(expm)
library(profvis)
library(coda)

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

lg_lklhd <- function(Q,input,output,time,etQ=NULL) {
  #if(! length(input)==length(output)) stop("Input and output must have same length.")
  S <- dim(Q)[1]
  #if(any(input>S)) stop("Input states cannot be larger than state space.")
  #if(any(output>S)) stop("Output states cannot be larger than state space.")
  N <- dim(input)[1]
  if(is.null(etQ)) {
    transMat <- expm::expm(time*Q)
  } else {
    transMat <- etQ
  }

  lklhd <- 0
  #transMat <- t(transMat)
  for(n in 1:N) {
  #   ntlStt <- rep(0,S)
  #   tptStt <- rep(0,S)
  #   ntlStt[input[n]] <- 1
  #   tptStt[output[n]] <- 1
   # lklhd <- lklhd + log(t(tptStt) %*% transMat %*% ntlStt)
  lklhd <- lklhd + log(input[n,] %*% transMat %*% output[,n])
  }
  # ntlStts <- matrix(0,N,S)
  # tptStts <- matrix(0,N,S)
  # ntlStts[cbind(1:N,input)] <- 1
  # tptStts[cbind(1:N,output)] <- 1
  # lklhd <- log(sum(diag(tptStts %*% t(transMat) %*% t(ntlStts))))
  
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
# ll <- lg_lklhd(Q,initialStates,t(output),1)
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

frstRdr <- function(t,Q,E,etQ){
  if(is.null(etQ)) {
    return(expm::expm(t*Q) %*% E * t)
  } else {
    return(etQ %*% E * t)
  }
}

frstRdrCrrctd <- function(t,Q,E) {
  out <- frstRdr(t,Q,E) - approx(t,Q,E)
  return(out)
}

lg_prr <- function(Q,tau,alpha) {
  priors <- - abs(log(Q)/tau)^alpha #dnorm(log(Q),log = TRUE)
  diag(priors) <- 0
  return(sum(priors)) # TODO: add bayesian bridge log prior
}

lg_prr_grd <- function(Q,tau,alpha) {
  # start with gaussian prior
  # gradient <- - log(Q)/tau^2
  # diag(gradient) <- 0
  
  # BB fixing a = 0.25
  lg_Q <- log(Q)
  gradient <- -  alpha * lg_Q / (tau^2) * abs(lg_Q/tau)^(alpha-2)
  diag(gradient) <- 0

  return(gradient) 
}

lg_lklhd_grd <- function(t,Q,input,output,method) {
  #if(! dim(input)[[==length(output)) stop("Input and output must have same length.")
  S <- dim(Q)[1]
  #if(any(input>S)) stop("Input states cannot be larger than state space.")
  #if(any(output>S)) stop("Output states cannot be larger than state space.")
  N <- dim(input)[[1]]
  if(! method %in% c("exact","approx","corrected")) stop("Use prescribed gradient method!")
  gradient <- Q * 0
  E <- Q * 0
  patty <- Q * 0
  
  outIn <- output %*% input
  
  for(i in 1:S) {
    E <- E*0
    E[i,i] <- 1
    # get diagonal contributions
    if (method=="exact") {
      diagContrib <- - blockwise(t,Q,E) * Q[i,i]
    } else if (method=="approx") {
      etQ <- expm::expm(Q*t)
      diagContrib <- - frstRdr(t,Q,E,etQ) * Q[i,i]
    } else {
      diagContrib <- - frstRdrCrrctd(t,Q,E) * Q[i,i]
    } 
    
    for(j in 1:S) {
      if(i != j) {
        patty <- patty * 0
        E <- E * 0
        E[i,j] <- 1
        if (method=="exact") {
          patty <- patty + blockwise(t,Q,E) * Q[i,j] + diagContrib
        } else if (method=="approx") {
          patty <- patty + frstRdr(t,Q,E,etQ) * Q[i,j] + diagContrib
        } else {
          patty <- patty + frstRdrCrrctd(t,Q,E) * Q[i,j] + diagContrib
        }
        #patty <- t(patty)
        #for(n in 1:N) {
        # ntlStt <- rep(0,S)
        # tptStt <- rep(0,S)
        # ntlStt[input[n]] <- 1
        # tptStt[output[n]] <- 1
        gradient[i,j] <- gradient[i,j] + sum(patty*outIn)   # gradient[i,j] + t(tptStt) %*% patty %*% ntlStt
        #}
      }
    }
  }
  
  # don't forget to divide by log-likelihood
  if (method=="approx") {
    gradient <- gradient / lg_lklhd(Q,input,output,t,etQ) 
  } else {
    gradient <- gradient / lg_lklhd(Q,input,output,t) 
  }
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
target <- function(Q,input,output,time,tau,alpha) {
  out <- lg_lklhd(Q,input,output,time) + lg_prr(Q,tau,alpha)
  return(out)
}

grad <- function(t,Q,input,output,method,tau,alpha) {
  out <- lg_lklhd_grd(t,Q,input,output,method) + lg_prr_grd(Q,tau,alpha)
  return(out)
}

delta <- function(n) {
  #return( min(0.01,n^(-0.5)) )
  return(n^(-0.1))
}

lg_Q_to_Q <- function(lg_Q) {
  Q <- exp(lg_Q)
  diag(Q) <- 0
  diag(Q) <- - rowSums(Q)
  return(Q)
}

hmc <- function(t,S,input,output,method,stepSize=0.01,L=20,maxIts=1000,
                targetAccept=0.8, Q_init=NULL, alpha=2) {
  
  #lg_Q <- list()
  if(is.null(Q_init)) {
    Q_init <- matrix(rexp(S^2),S,S)
    diag(Q_init) <- 0
    diag(Q_init) <- - rowSums(Q_init) 
  }
  lg_Q <- log(Q_init)
  diag(lg_Q) <- 0
  cat(lg_Q, "\n", file="data/highDimSim_20_half.txt",append=TRUE)
  
  #diag(lg_Q) <- 0
  
  taus <- rep(0,maxIts) # global parameters on bayesian bridge prior
  taus[1] <- 1 #0.0001 # tau^(-0.25) ~ Gamma(shape=1,rate=2) prior
  
  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  
  currentU  <- - target(Q_init,input,output,t,taus[1],alpha)
  log_probs <- rep(0,maxIts)
  log_probs[1] <- currentU
  
  for (i in 2:maxIts) {
    proposalState    <- lg_Q
    momentum         <- matrix(rnorm(S^2),S,S)
    diag(momentum)   <- 0
    currentK         <- sum(momentum^2)/2
    
    # leapfrog steps
    momentum <- momentum + 0.5 * stepSize * grad(t,lg_Q_to_Q(proposalState),input,output,method,taus[i-1],alpha)
    for (l in 1:L) {
      proposalState <- proposalState + stepSize * momentum
      # diag(proposalState) <- 0
      # diag(proposalState) <- - rowSums(proposalState) 
      if (l!=L) momentum <- momentum + stepSize * grad(t,lg_Q_to_Q(proposalState),input,output,method,taus[i-1],alpha)
    }
    momentum <- momentum + 0.5 * stepSize * grad(t,lg_Q_to_Q(proposalState),input,output,method,taus[i-1],alpha)
    
    # quantities for accept/reject
    proposedU = - target(lg_Q_to_Q(proposalState),input,output,t,taus[i-1],alpha)
    proposedK = sum(momentum^2)/2
    u <- runif(1)
    
    if (log(u) < currentU + currentK - proposedU - proposedK) {
      lg_Q    <- proposalState
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
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
    
    # nu <- rgamma(n=1,
    #              shape=1+(S*(S-1))*4,
    #              rate=2+sum(abs(lg_Q[[i]])^(1/4)))
    #taus[i] <- nu^(-4)
    taus[i] <-rgamma(n=1,
                      shape=1+(S*(S-1))/alpha,
                      rate=2+sum(abs(lg_Q)^(alpha),na.rm = TRUE))
    taus[i] <- taus[i]^(-1/alpha)
    currentU <- - target(lg_Q_to_Q(lg_Q),input,output,t,taus[i],alpha)
    
    if (i %% 100 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n") 
    
    if(i %% 100 == 0) cat(lg_Q, "\n", file="data/highDimSim_20_half.txt",append=TRUE)
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(list(log_probs=-log_probs,taus=taus))
}

sim_bridge <- function(n=100,tau=1,alpha=1) {
  # gibbs sampler for sampling bridge variates
  omegas <- (1+alpha)/2*rgamma(n=n,shape=2+1/alpha,rate=1) + (1-alpha)/2*rgamma(n=n,shape=1+1/alpha,rate=1)
  betas  <- runif(n=n,min=-0.1,max=0.1) #rnorm(n=n,sd=1)
  us     <- runif(n=n, max=1-abs(betas/tau)*omegas^(-1/alpha))
  
  #bs     <- 1/tau * (1-us) * omegas ^
  for(i in 1:1000) {
    omegas <- (1+alpha)/2*rgamma(n=n,shape=2+1/alpha,rate=1) + (1-alpha)/2*rgamma(n=n,shape=1+1/alpha,rate=1)
    bs     <-  1/tau * (1-us) * omegas^(1/alpha)
    betas  <- runif(n=n,min=-bs,max=bs)
    us     <- runif(n=n, max=1-abs(betas/tau)*omegas^(-1/alpha))
  }
  return(betas)
}

# test
set.seed(1)
S <- 20
maxIts <- 10000000
N <- 300
entries <- sim_bridge(n=S^2,tau=1,alpha=0.5) 
entries <- entries / max(abs(entries))
initialStates <- sample(x=1:S,size=N,replace=TRUE)
Q       <- matrix(entries,S,S)
Q       <- exp(Q)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
output <- simCTMC(Q,initialStates,1)

ntlStts <- matrix(0,N,S)
tptStts <- matrix(0,N,S)
ntlStts[cbind(1:N,initialStates)] <- 1
tptStts[cbind(1:N,output)] <- 1
tptStts <- t(tptStts)

#profvis({
hmcOut <- hmc(t=1,S=S,input=ntlStts,output=tptStts,stepSize=0.001,
              method="approx",L=4,maxIts = maxIts, targetAccept = 0.65, Q_init=Q,
              alpha=0.5)
#})
#saveRDS(hmcOut,file="data/lowDimSimApprox.rds")





# 
# plot(hmcOut[[2]],type="l")
# effectiveSize(hmcOut[[2]]) # 11865.41
# 
# chain <- hmcOut[[1]]
# 
# variable <- rep(0,maxIts)
# for(i in 1:maxIts) {
#   variable[i] <- chain[[i]][2,1]
# }
# plot(variable,type="l")
# abline(h=log(Q[2,1]),col="red")
# effectiveSize(variable) #231.2789
# 
# variable <- rep(0,maxIts)
# for(i in 1:maxIts) {
#   variable[i] <- chain[[i]][5,3]
# }
# plot(variable,type="l")
# abline(h=log(Q[5,3]),col="red")
# effectiveSize(variable) #301
# 
# plot(hmcOut[[3]][2:maxIts],type="l")
# abline(h=1*sqrt(2),col="red")
# effectiveSize(hmcOut[[3]][2:maxIts]) # 79

