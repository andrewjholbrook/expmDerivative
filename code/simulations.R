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
  transMat <- expm(time*Q)
  
  lklhd <- 0
  for(n in 1:N) {
    ntlStt <- rep(0,S)
    tptStt <- rep(0,S)
    ntlStt[input[n]] <- 1
    tptStt[output[n]] <- 1
    lklhd <- lklhd + log(t(tptStt) %*% t(transMat) %*% ntlStt)
  }
  return(lklhd)
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



