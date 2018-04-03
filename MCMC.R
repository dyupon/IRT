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
    cond.prob.current <- (colSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- (colSums(log(P.proposal) * (1 - X) + log(1 - P.proposal)*(1 - X)))
    ratio <- exp(cond.prob.proposal + log(theta.prior(val = proposal)) -
                        (cond.prob.current + log(theta.prior(val = chain.theta))))
    stopifnot(!any(is.na(ratio)))
    
  } else {
    n <- nitems
    P.proposal <- outer(proposal, chain.theta, probability)
    cond.prob.current <- exp(rowSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(rowSums(log(P.proposal) * X + log(1 - P.proposal)*(1 - X)))
    ratio <- cond.prob.proposal*b.prior(val = proposal) /
                        (cond.prob.current*b.prior(val = chain.beta))
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
iterations <- 10000

mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
fitdistr(mc$chain.theta[iterations,], "normal")
fitdistr(mc$chain.beta[iterations,], "normal")
hist(mc$chain.theta[, 1],nclass = 30, main = "Posterior of theta_1")
abline(v = mean(mc$chain.theta[,1]), col = "red")
hist(mc$chain.beta[,1],nclass = 30, main = "Posterior of beta_1")
abline(v = mean(mc$chain.beta[,1]), col = "red")
plot(as.vector(mc$chain.theta[,1]), type = "l", main = "Chain values of theta_1" )
plot(as.vector(mc$chain.beta[,1]), type = "l",  main = "Chain values of beta_1" )
hist(mc$chain.theta[, 50],nclass = 30, main = "Posterior of theta_50")
abline(v = mean(mc$chain.theta[,50]), col = "red")
plot(as.vector(mc$chain.theta[,50]), type = "l", main = "Chain values of theta_50" )