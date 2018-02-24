library(ggplot2)
library(MCMCpack)
set.seed(62)
#DATA SIMULATION
rasch.modeling <- function(persons, items) {
  temp <- matrix(rep(items, length(persons)), ncol = length(persons))
  exp_index <- t(apply(temp, 1, '-', persons))
  probabilities <- 1 / (1 + exp(exp_index))
  output <- list()
  output$items <- items
  output$persons <- persons
  output$data <- matrix(apply(probabilities, 1, function(x) as.numeric(runif(1) > x)), ncol = length(persons))
  output
}
npersons <- 100
nitems <- 20
items <- rnorm(nitems, 0, 1)
persons <- rnorm(npersons)
rmdata <- rasch.modeling(persons, items)
X <- rmdata$data

probability <- function(beta, theta) {
  return(1 / (1 + exp(beta - theta)))
}

# loglikelihood of the data
data.like <- function(param, data){
  theta <- param[1] # person ability
  b <- param[2] # item difficulty
  p <- exp(theta - b)/(1 + exp(theta - b))
  y <- data # data
  # log likelihood
  like <- sum(y*log(p) + (1 - y)*log(1 - p)) 
  return(like)
}

theta.prior <- function(mean = 0, sd = 1, val){
  mu <- mean # prior mean
  s <- sd # prior standard deviation
  y <- val # observed "x" value from prior
  prior <- dnorm(x = y, mean = mu, sd = s) 
  return(prior)
}

b.prior <- function(mean = 0, sd = 1, val){
  mu <- mean
  s <- sd
  y <- val
  prior <- dnorm(x = y, mean = mu, sd = s)
  return(prior)
} 

# proposal, centered around the current value of the chain
prop.theta <- function(current){
  proposal <- rnorm(current, mean = current, sd = 1)
  return(proposal)
}

prop.beta <- function(current){
  proposal <- rnorm(current, mean = current, sd = 1)
  return(proposal)
}

acceptance <- function(proposal, chain.theta, chain.beta, X, mode) {
  P <- outer(chain.beta, chain.theta, probability)
  nitems <- length(chain.beta)
  npersons <- length(chain.theta)
  if (mode == 0) {
    n <- npersons
    P.proposal <- outer(chain.beta, proposal, probability)
    cond.prob.current <- exp(colSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(colSums(log(P.proposal) * X + log(1 - P.proposal)*(1 - X)))
  } else {
    n <- nitems
    P.proposal <- outer(proposal, chain.theta, probability)
    cond.prob.current <- exp(rowSums(log(P) * X + log(1 - P)*(1 - X)))
    cond.prob.proposal <- exp(rowSums(log(P.proposal) * X + log(1 - P.proposal)*(1 - X)))
  }
  ratio <- cond.prob.proposal*theta.prior(val = proposal)/(cond.prob.current*b.prior(val = chain.beta))
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
    chain.theta[i,] <- ifelse(acceptance(t.proposal, chain.theta[i - 1,], chain.beta[i - 1,], X, 0) < runif(1), 
           t.proposal,
           chain.theta[i - 1,])
    chain.beta[i,] <- ifelse(acceptance(b.proposal, chain.theta[i - 1,], chain.beta[i - 1,], X, 1) < runif(1),
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

par(mfrow = c(1,1))
burnIn <- 25
hist(mc$chain.theta[, 1],nclass = 30, main = "Posterior of theta_1")
abline(v = mean(mc$chain.theta[,1]), col = "red")
hist(mc$chain.beta[,1],nclass = 30, main = "Posterior of beta_1")
abline(v = mean(mc$chain.beta[,1]), col = "red")
plot(as.vector(mc$chain.theta[,1]), type = "l", main = "Chain values of theta_1" )
plot(as.vector(mc$chain.beta[,1]), type = "l",  main = "Chain values of beta_1" )
print(mc$chain.theta[1:burnIn, 1])
print(mc$chain.beta[1:burnIn, 1])
fitdistr(mc$chain.theta[burnIn,], "normal")
fitdistr(mc$chain.beta[burnIn,], "normal")
