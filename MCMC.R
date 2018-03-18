library(ggplot2)
library(MCMCpack)
set.seed(62)
#DATA SIMULATION
sim.rasch <- function(persons, items, cutpoint = "randomized"){
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
  #print("Next")
  #print("current state")
  #print(current)
  #print("proposal state")
  #print(proposal)
  stopifnot(is.finite(proposal))
  return(proposal)
}

acceptance <- function(proposal, chain.theta, chain.beta, X, mode) {
  P <- outer(chain.beta, chain.theta, probability)
  if (mode == 0) {
    n <- npersons
    P.proposal <- outer(chain.beta, proposal, probability)
    cond.prob.current <- exp(colSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(colSums(log(P.proposal) * X + log(1 - P.proposal)*(1 - X)))
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
  X <- data
  chain.beta <- array(dim = c(iterations + 1, nitems))
  chain.theta <- array(dim = c(iterations + 1, npersons))
  chain.beta[1,] <- startvalue.beta 
  chain.theta[1,] <- startvalue.theta
  # iterations
  for (i in 2:iterations) {
    t.proposal <- prop.theta(chain.theta[i - 1,]) # sample a new value for theta
    b.proposal <- prop.beta(chain.beta[i - 1,])
    #print("chain beta i-1")
    #print(chain.beta[i - 1,])
    acceptance.theta <- acceptance(t.proposal, chain.theta[i - 1,], chain.beta[i - 1,], X, 0)
    acceptance.theta[is.na(acceptance.theta)] <- 0
    #print("acceptance theta")
    #print(acceptance.theta)
    chain.theta[i,] <- ifelse(acceptance.theta == rep(1, length(acceptance.theta)), 
                              t.proposal,
                              chain.theta[i - 1,])
    chain.theta[i,] <- ifelse(acceptance.theta < runif(length(acceptance.theta), 0, 1), 
                              t.proposal,
                              chain.theta[i - 1,])
    acceptance.beta <- acceptance(b.proposal, chain.theta[i,], chain.beta[i - 1,], X, 1)
    acceptance.beta[is.na(acceptance.beta)] <- 0
    #print("acceptance beta")
    #print(acceptance.beta)
    chain.beta[i,] <- ifelse(acceptance.beta == rep(1, length(acceptance.beta)),
                             b.proposal,
                             chain.beta[i - 1,])
    chain.beta[i,] <- ifelse(acceptance.beta < runif(length(acceptance.beta), 0, 1),
                             b.proposal,
                             chain.beta[i - 1,])
    #print(b.proposal)
    #print("chain i")
    #print(chain.beta[i,])
    #stopifnot(is.finite(t.proposal))
      }
  list <- list("chain.theta" = mcmc(chain.theta), "chain.beta" = mcmc(chain.beta))
  return(list)
}

startvalue.theta <- rep(0.1, npersons)
startvalue.beta <- rep(0.1, nitems)
iterations <- 50000

mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
burnIn <- 2500
hist(mc$chain.theta[burnIn:iterations - 1, 1],nclass = 30, main = "Posterior of theta_1")
abline(v = mean(mc$chain.theta[burnIn:iterations - 1,1]), col = "red")
hist(mc$chain.beta[,1],nclass = 30, main = "Posterior of beta_1")
abline(v = mean(mc$chain.beta[burnIn:iterations - 1,1]), col = "red")
plot(as.vector(mc$chain.theta[burnIn:iterations - 1,1]), type = "l", main = "Chain values of theta_1" )
plot(as.vector(mc$chain.beta[burnIn:iterations - 1,1]), type = "l",  main = "Chain values of beta_1" )
fitdistr(mc$chain.theta[burnIn:iterations - 1,], "normal")
fitdistr(mc$chain.beta[burnIn:iterations - 1,], "normal")
autocorr(mcmc(mc$chain.beta[1:iterations - 1,1]), lags <- c(0, 10000, 20000, 30000, 40000, 49000))
autocorr(mcmc(mc$chain.theta[1:iterations - 1,1]), lags <- c(0, 10000, 20000, 30000, 40000, 49000))
