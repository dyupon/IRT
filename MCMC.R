library(ggplot2)
library(MCMCpack)
library(matrixStats)
set.seed(62)

sim.rasch <- function(persons, items){
  schwierig <- items
  n.items <- length(items)
  faehig <- persons
  n.persons <- length(persons)
  fsmat <- outer(faehig, schwierig, "-")
  psolve <- exp(fsmat)/(1 + exp(fsmat))
  R <- (matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
  R
}

npersons <- 100
nitems <- 20
items <- rnorm(nitems, 0, 1)
persons <- rnorm(npersons)
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)

probability <- function(beta, theta) {
  return(1 / (1 + exp(beta - theta)))
}

theta.prior <- function(mean = 0, sd = 1, val){
  mu <- mean 
  s <- sd 
  y <- val
  prior <- dnorm(x = y, mean = mu, sd = s) 
  prior
}

b.prior <- function(mean = 0, sd = 1, val){
  mu <- mean
  s <- sd
  y <- val
  prior <- dnorm(x = y, mean = mu, sd = s)
  prior
} 

# proposal, centered around the current value of the chain
prop.theta <- function(current){
  proposal <- rnorm(current, mean = current, sd = 1)
  proposal
}

prop.beta <- function(current){
  proposal <- rnorm(current, mean = current, sd = 1)
  return(proposal)
}

acceptance <- function(proposal, chain.theta, chain.beta, X, mode) {
  P <- outer(chain.beta, chain.theta, probability)
  if (mode == 0) {
    n <- npersons
    P.proposal <- outer(chain.beta, proposal, probability)
    cond.prob.current <- exp(colSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(colSums(log1p(P.proposal) * (1 - X) + log(1 - P.proposal)*(1 - X)))
    ratio <- cond.prob.proposal*theta.prior(val = proposal)/(cond.prob.current*b.prior(val = chain.beta))
  } else {
    n <- nitems
    P.proposal <- outer(proposal, chain.theta, probability)
    cond.prob.current <- exp(rowSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(rowSums(log(P.proposal) * X + log(1 - P.proposal)*(1 - X)))
    ratio <- cond.prob.proposal*b.prior(val = proposal)/(cond.prob.current*b.prior(val = chain.beta))
  }
  return(pmin(ratio, rep(1,n)))
}


#GS alg
rasch.gs <- function(startvalue.theta, startvalue.beta, iterations, data, nitems, npersons) {
  chain.beta <- array(dim = c(iterations + 1, nitems))
  chain.theta <- array(dim = c(iterations + 1, npersons))
  chain.beta[1,] <- startvalue.beta 
  chain.theta[1,] <- startvalue.theta
  # iterations
  for (i in 2:(iterations + 1)) {
    t.proposal <- prop.theta(chain.theta[i - 1,])
    b.proposal <- prop.beta(chain.beta[i - 1,])
    acceptance.theta <- acceptance(t.proposal, chain.theta[i - 1,], 
                                   chain.beta[i - 1,], X, 
                                   0)
    chain.theta[i,] <- ifelse(runif(length(acceptance.theta), 0, 1) < acceptance.theta, 
                              t.proposal,
                              chain.theta[i - 1,])
    acceptance.beta <- acceptance(b.proposal, chain.theta[i,], 
                                  chain.beta[i - 1,], X, 
                                  1)
    chain.beta[i,] <- ifelse(runif(length(acceptance.beta), 0, 1) < acceptance.beta,
                             b.proposal,
                             chain.beta[i - 1,])
  }
  list <- list("chain.theta" = mcmc(chain.theta), "chain.beta" = mcmc(chain.beta))
  return(list)
}

startvalue.theta <- rep(0.1, npersons)
startvalue.beta <- rep(0.1, nitems)
iterations <- 5000

mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
