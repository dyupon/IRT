library(ggplot2)
library(MCMCpack)
library(matrixStats)
library(beepr)
library(tidyr)
library(dplyr)
library(latex2exp)
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


probability <- function(beta, theta) {
  1 / (1 + exp(beta - theta))
}

theta.prior <- function(mean = 0, sd = 1, val){
  dnorm(x = val, mean, sd, log = T) 
}

b.prior <- function(mean = 0, sd = 1, val){
  dnorm(x = val, mean, sd, log = T)
} 

prop.theta <- function(current){
  current + rnorm(n = current, 0, sd = 1)
}

prop.beta <- function(current){
  current + rnorm(n = current, 0, sd = 1)
}

acceptance <- function(proposal, chain.theta, chain.beta, X, mode) {
  P <- outer(chain.beta, chain.theta, probability)
  diff <- outer(chain.beta, chain.theta, function(x, y) {x - y})
  lgP <- log(P)
  if (mode == 0) {
    n <- npersons
    P.proposal <- outer(chain.beta, proposal, probability)
    diff.proposal <- outer(chain.beta, proposal, function(x, y) {x - y})
    lgP.proposal <- log(P.proposal)
    cond.prob.current <- colSums(lgP * X + (diff + lgP)*(1 - X))
    cond.prob.proposal <- colSums(lgP.proposal * X + (diff.proposal + lgP.proposal)*(1 - X))
    ratio <- cond.prob.proposal + theta.prior(val = proposal) - 
      (cond.prob.current + theta.prior(val = chain.theta))
    stopifnot(!any(is.na(ratio)))
  } else {
    n <- nitems
    P.proposal <- outer(proposal, chain.theta, probability)
    diff.proposal <- outer(proposal, chain.theta, function(x, y) {x - y})
    lgP.proposal <- log(P.proposal)
    cond.prob.current <- rowSums(lgP * X + (diff + lgP)*(1 - X))
    cond.prob.proposal <- rowSums(lgP.proposal * X + (diff.proposal + lgP.proposal)*(1 - X))
    ratio <- cond.prob.proposal + b.prior(val = proposal) -
      (cond.prob.current + b.prior(val = chain.beta))
  }
  return(pmin(ratio, rep(0,n)))
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
    chain.theta[i,] <- ifelse(-rexp(n = length(acceptance.theta), 1) < acceptance.theta, 
                              t.proposal,
                              chain.theta[i - 1,])
    acceptance.beta <- acceptance(b.proposal, chain.theta[i,], 
                                  chain.beta[i - 1,], X, 
                                  1)
    chain.beta[i,] <- ifelse(-rexp(n = length(acceptance.beta), 1) < acceptance.beta,
                             b.proposal,
                             chain.beta[i - 1,])
  }
  list <- list("chain.theta" = mcmc(chain.theta), "chain.beta" = mcmc(chain.beta))
  return(list)
}

get.graphs <- function(mc, npersons, items, iterations, trueMeanPersons, trueMeanItems) {
  df.theta <- as.data.frame(mc$chain.theta)
  colnames(df.theta)[1] <- "theta1"
  colnames(df.theta)[npersons] <- paste0("theta", npersons)
  last.state.theta <- as.data.frame(t(df.theta[iterations + 1,]))
  colnames(last.state.theta) <- "last.state"
  df.theta %>% gather("theta", "value", c(theta1, paste0("theta", npersons))) %>% select(theta, value) -> posterior
  posterior %>% group_by(theta) %>% summarize(group.mean = mean(value)) -> mean.posterior
  chain.val.theta <- cbind.data.frame(mc$chain.theta[,1], mc$chain.theta[,(npersons/2)], mc$chain.theta[,npersons])
  colnames(chain.val.theta) <- c("theta1", paste0("theta", npersons/2), paste0("theta", npersons))
  chain.val.theta %>% mutate(index = as.numeric(row.names(chain.val.theta))) -> chain.val.theta
  df.beta <- as.data.frame(mc$chain.beta)
  colnames(df.beta)[1] <- "beta1"
  colnames(df.beta)[nitems] <- paste0("beta", nitems)
  last.state.beta <- as.data.frame(t(df.beta[iterations + 1,]))
  colnames(last.state.beta) <- "last.state"
  df.beta %>% gather("beta", "value", c(beta1, paste0("beta", nitems))) %>% select(beta, value) -> posterior.beta
  posterior.beta %>% group_by(beta) %>% summarize(group.mean = mean(value)) -> mean.posterior.beta
  chain.val.beta <- cbind.data.frame(mc$chain.beta[,1], mc$chain.beta[,(nitems/2)], mc$chain.beta[,nitems])
  colnames(chain.val.beta) <- c("beta1", paste0("beta", nitems/2), paste0("beta", nitems))
  chain.val.beta %>% mutate(index = as.numeric(row.names(chain.val.beta))) -> chain.val.beta
  
  p1 <- ggplot(data = last.state.theta, aes(x = last.state)) +
    geom_histogram(aes(x = last.state), binwidth = 0.1, color = "red",position = "identity", alpha = 0.1) +
    geom_vline(aes(xintercept = mean(last.state)),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = TeX("Posterior of $\\theta = (\\theta_{1}, ..., \\theta_{N})$ at the last state of chain"), y = "Frequency", x = "Value") +
    geom_vline(aes(xintercept = trueMeanPersons),
               color = "green", linetype = "dashed", size = 1) 
  
  p2 <- ggplot(data = posterior, aes(x = value, fill = theta, color = theta)) +
    geom_histogram(aes(x = value), binwidth = 0.1,position = "identity", alpha = 0.1) +
    geom_vline(x = mean.posterior, aes(xintercept = mean.posterior[1,2]),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(x = mean.posterior, aes(xintercept = mean.posterior[2,2]),
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = TeX("Trajectories of $\\theta_{1}$ and $\\theta_{1000}$"),y = "Frequency", x = "Value") +
    theme(legend.title = element_blank())
  
  p3 <- ggplot(data = chain.val.theta, aes(x = index, y = theta1)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{1}$"), x = "The number of chain state", y = "Value")
  p4 <- ggplot(data = chain.val.theta, aes(x = index, y = theta500)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{500}$"), x = "The number of chain state", y = "Value")
  p5 <- ggplot(data = chain.val.theta, aes(x = index, y = theta1000)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{1000}$"), x = "The number of chain state", y = "Value")
  
  
  p6 <- ggplot(data = last.state.beta, aes(x = last.state)) +
    geom_histogram(aes(x = last.state), binwidth = 0.1, color = "red",position = "identity", alpha = 0.1) +
    geom_vline(aes(xintercept = mean(last.state)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = trueMeanItems),
               color = "green", linetype = "dashed", size = 1) +
    labs(title = TeX("Posterior of $\\beta = (\\beta_{1}, ..., \\beta_{J})$ at the last state of chain"), y = "Frequency", x = "Last state")
  
  p7 <- ggplot(data = posterior.beta, aes(x = value, fill = beta, color = beta)) +
    geom_histogram(aes(x = value), binwidth = 0.1,position = "identity", alpha = 0.1) +
    geom_vline(x = mean.posterior.beta, aes(xintercept = mean.posterior.beta[1,2]),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(x = mean.posterior.beta, aes(xintercept = mean.posterior.beta[2,2]),
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = TeX("Trajectories of $\\beta_{1}$ and $\\beta_{40}$"), y = "Frequency", x = "Value") +
    theme(legend.title = element_blank())
  
  p8 <- ggplot(data = chain.val.beta, aes(x = index, y = beta1)) + geom_line() +
    labs(title = TeX("Chain values of $\\beta_{1}$"), x = "The number of chain state", y = "Value")
  p9 <- ggplot(data = chain.val.beta, aes(x = index, y = beta20)) + geom_line() +
    labs(title = TeX("Chain values of $\\beta_{20}$"), x = "The number of chain state", y = "Value")
  p10 <- ggplot(title = TeX("Chain values of $\\beta_{40}$"), data = chain.val.beta, aes(x = index, y = beta40)) + geom_line() +
    labs(title = TeX("Chain values of $\\beta_{40}$"), x = "The number of chain state", y = "Value")
  
  lags <- seq(0, 100, 1)
  autocorr.theta <- cbind.data.frame(autocorr(mc$chain.theta[,1], lags = lags), 
                                     autocorr(mc$chain.theta[,npersons/2], lags = lags),
                                     autocorr(mc$chain.theta[,npersons], lags = lags))
  colnames(autocorr.theta) <- c("theta1", paste0("theta", npersons/2), paste0("theta", npersons))
  autocorr.theta %>% mutate(index = lags) -> autocorr.theta
  autocorr.theta  %>% gather("theta", "value", 1:(dim(autocorr.theta)[2] - 1)) -> autocorr.theta
  p11 <- ggplot(data = autocorr.theta, aes(x = index, y = value, color = theta)) + 
    geom_line() + 
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    theme(legend.title = element_blank()) +
    labs(title = TeX("Autocorrelation function for chain $\\theta_{1}, \\theta_{500}, \\theta_{1000}$"), y = "Autocorrelation", x = "Lags")
  
  autocorr.beta <- cbind.data.frame(autocorr(mc$chain.beta[,1], lags = lags), 
                                    autocorr(mc$chain.beta[,nitems/2], lags = lags),
                                    autocorr(mc$chain.beta[,nitems], lags = lags))
  colnames(autocorr.beta) <- c("beta1", paste0("beta", nitems/2), paste0("beta", nitems))
  autocorr.beta %>% mutate(index = lags) -> autocorr.beta
  autocorr.beta  %>% gather("beta", "value", 1:(dim(autocorr.beta)[2] - 1)) -> autocorr.beta
  p12 <- ggplot(data = autocorr.beta, aes(x = index, y = value, color = beta)) + 
    geom_line() + 
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    theme(legend.title = element_blank()) +
    labs(title = TeX("Autocorrelation function for chain $\\beta_{1}, \\beta_{20}, \\beta_{40}$"), y = "Autocorrelation", x = "Lags")
  return(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12))
}

get.graphs1 <- function(mc1, npersons, nitems, iterations) {
  df.theta.1 <- as.data.frame(mc1$chain.theta)
  colnames(df.theta.1)[1] <- "theta1"
  colnames(df.theta.1)[npersons] <- paste0("theta", npersons)
  last.state.theta.1 <- as.data.frame(t(df.theta.1[iterations + 1,]))
  colnames(last.state.theta.1) <- "last.state"
  df.theta.1 %>% gather("theta", "value", c(theta1, paste0("theta", npersons))) %>% select(theta, value) -> posterior.1
  
  posterior.1 %>% group_by(theta) %>% summarize(group.mean = mean(value)) -> mean.posterior.1
  chain.val.theta.1 <- cbind.data.frame(mc1$chain.theta[,1], mc1$chain.theta[,(npersons/2)], mc1$chain.theta[,npersons])
  colnames(chain.val.theta.1) <- c("theta1", paste0("theta", npersons/2), paste0("theta", npersons))
  chain.val.theta.1 %>% mutate(index = as.numeric(row.names(chain.val.theta.1))) -> chain.val.theta.1
  
  p1 <- ggplot(data = last.state.theta.1, aes(x = last.state)) +
    geom_histogram(aes(x = last.state), binwidth = 0.1, color = "red",position = "identity", alpha = 0.1) +
    labs(title = TeX("Posterior of $\\theta = (\\theta_{1}, ..., \\theta_{N})$ at the last state of chain"), y = "Frequency", x = "Value") +
    geom_vline(aes(xintercept = trueMeanPersons),
               color = "green", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = -trueMeanPersons),
               color = "green", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(mc1$chain.theta[iterations,1:npersons/2])),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(mc1$chain.theta[iterations,(npersons/2):npersons])),
               color = "red", linetype = "dashed", size = 1) 
  
  p2 <- ggplot(data = posterior.1, aes(x = value, fill = theta, color = theta)) +
    geom_histogram(aes(x = value), binwidth = 0.1,position = "identity", alpha = 0.1) +
    geom_vline(x = mean.posterior.1, aes(xintercept = mean.posterior.1[1,2]),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(x = mean.posterior.1, aes(xintercept = mean.posterior.1[2,2]),
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = TeX("Trajectories of $\\theta_{1}$ and $\\theta_{1000}$"),y = "Frequency", x = "Value") +
    theme(legend.title = element_blank())
  
  p3 <- ggplot(data = chain.val.theta.1, aes(x = index, y = theta1)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{1}$"), x = "The number of chain state", y = "Value")
  p4 <- ggplot(data = chain.val.theta.1, aes(x = index, y = theta500)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{500}$"), x = "The number of chain state", y = "Value")
  p5 <- ggplot(data = chain.val.theta.1, aes(x = index, y = theta1000)) + geom_line() +
    labs(title = TeX("Chain values of $\\theta_{1000}$"), x = "The number of chain state", y = "Value")
  
  lags <- seq(0, 100, 1)
  autocorr.theta.1 <- cbind.data.frame(autocorr(mc1$chain.theta[,1], lags = lags), 
                                       autocorr(mc1$chain.theta[,npersons/2], lags = lags),
                                       autocorr(mc1$chain.theta[,npersons], lags = lags))
  colnames(autocorr.theta.1) <- c("theta1", paste0("theta", npersons/2), paste0("theta", npersons))
  autocorr.theta.1 %>% mutate(index = lags) -> autocorr.theta.1
  autocorr.theta.1  %>% gather("theta", "value", 1:(dim(autocorr.theta.1)[2] - 1)) -> autocorr.theta.1
  
  p6 <- ggplot(data = autocorr.theta.1, aes(x = index, y = value, color = theta)) + 
    geom_line() + 
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    theme(legend.title = element_blank()) +
    labs(title = TeX("Autocorrelation function for chain $\\theta_{1}, \\theta_{500}, \\theta_{1000}$"), y = "Autocorrelation", x = "Lags")
  return(list(p1, p2, p3, p4, p5, p6))
}

print.summary <- function(mc, iterations, persons) {
  print("fitdistr theta_1 ... theta_N at the last chain state to normal")
  print(fitdistr(mc$chain.theta[iterations,], "normal"))
  print(min(mc$chain.theta[iterations,]))
  print(max(mc$chain.theta[iterations,]))
  print(min(persons))
  print(max(persons))
  print("fitdistr beta_1 ... beta_J at the last chain state to normal")
  print(fitdistr(mc$chain.beta[iterations,], "normal"))
}

print.summary.1 <- function(mc1, iterations, npersons) {
  print(fitdistr(mc1$chain.theta[iterations,1:npersons/2], "normal"))
  print(fitdistr(mc1$chain.theta[iterations,(npersons/2):npersons], "normal"))
}
npersons <- 1000
nitems <- 40
startvalue.theta <- rep(0.1, npersons)
startvalue.beta <- rep(0.1, nitems)
iterations <- 10000

trueMeanPersons <- 5

######## N(0, 1) ######## 
trueMeanItems <- 0
persons <- rnorm(npersons,trueMeanPersons, 1)
items <- rnorm(nitems, trueMeanItems, 1)
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)
mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
get.graphs(mc, npersons, items, iterations, trueMeanPersons, trueMeanItems)
print.summary(mc, iterations, persons)

######## N(0, 5) ######## 
iterations <- 100000
persons <- rnorm(npersons,trueMeanPersons, 1)
items <- rnorm(nitems, trueMeanItems, 5)
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)
mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
get.graphs(mc, npersons, items, iterations, trueMeanPersons, trueMeanItems)
print.summary(mc, iterations, persons)

######## N(3, 1) ######## 
iterations <- 10000
trueMeanItems <- 3
persons <- rnorm(npersons,trueMeanPersons, 1)
items <- rnorm(nitems, trueMeanItems, 1)
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)
mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
get.graphs(mc, npersons, items, iterations, trueMeanPersons, trueMeanItems)
print.summary(mc, iterations, persons)

######## N(-3, 1) ######## 
iterations <- 10000
trueMeanItems <- -3
persons <- rnorm(npersons,trueMeanPersons, 1)
items <- rnorm(nitems, trueMeanItems, 1)
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)
mc <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
get.graphs(mc, npersons, items, iterations, trueMeanPersons, trueMeanItems)
print.summary(mc, iterations, persons)


######## MEAN COMPARISON ######## 
iterations <- 20000
trueMeanItems <- 2
items <- rnorm(nitems, 2, 4)
persons <- c(rnorm(npersons/2, -trueMeanPersons, 1), rnorm(npersons/2, trueMeanPersons, 1))
rmdata <- sim.rasch(persons, items)
X <- t(rmdata)
mc1 <- rasch.gs(startvalue.theta, startvalue.beta, iterations, X, nitems, npersons)
get.graphs1(mc1, npersons, nitems, iterations)
print.summary.1(mc1, iterations, npersons)
df.beta <- as.data.frame(mc1$chain.beta)
colnames(df.beta)[1] <- "beta1"
colnames(df.beta)[nitems] <- paste0("beta", nitems)
last.state.beta <- as.data.frame(t(df.beta[iterations + 1,]))
colnames(last.state.beta) <- "last.state"
df.beta %>% gather("beta", "value", c(beta1, paste0("beta", nitems))) %>% select(beta, value) -> posterior.beta
posterior.beta %>% group_by(beta) %>% summarize(group.mean = mean(value)) -> mean.posterior.beta
chain.val.beta <- cbind.data.frame(mc1$chain.beta[,1], mc1$chain.beta[,(nitems/2)], mc1$chain.beta[,nitems])
colnames(chain.val.beta) <- c("beta1", paste0("beta", nitems/2), paste0("beta", nitems))
chain.val.beta %>% mutate(index = as.numeric(row.names(chain.val.beta))) -> chain.val.beta
p6 <- ggplot(data = last.state.beta, aes(x = last.state)) +
  geom_histogram(aes(x = last.state), color = "red",position = "identity", alpha = 0.1) +
  geom_vline(aes(xintercept = mean(last.state)),
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = trueMeanItems),
             color = "green", linetype = "dashed", size = 1) +
  labs(title = TeX("Posterior of $\\beta = (\\beta_{1}, ..., \\beta_{J})$ at the last state of chain"), y = "Frequency", x = "Last state")

p7 <- ggplot(data = posterior.beta, aes(x = value, fill = beta, color = beta)) +
  geom_histogram(aes(x = value),position = "identity", alpha = 0.1) +
  geom_vline(x = mean.posterior.beta, aes(xintercept = mean.posterior.beta[1,2]),
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(x = mean.posterior.beta, aes(xintercept = mean.posterior.beta[2,2]),
             color = "blue", linetype = "dashed", size = 1) +
  labs(title = TeX("Trajectories of $\\beta_{1}$ and $\\beta_{40}$"), y = "Frequency", x = "Value") +
  theme(legend.title = element_blank())

p8 <- ggplot(data = chain.val.beta, aes(x = index, y = beta1)) + geom_line() +
  labs(title = TeX("Chain values of $\\beta_{1}$"), x = "The number of chain state", y = "Value")
p9 <- ggplot(data = chain.val.beta, aes(x = index, y = beta20)) + geom_line() +
  labs(title = TeX("Chain values of $\\beta_{20}$"), x = "The number of chain state", y = "Value")
p10 <- ggplot(title = TeX("Chain values of $\\beta_{40}$"), data = chain.val.beta, aes(x = index, y = beta40)) + geom_line() +
  labs(title = TeX("Chain values of $\\beta_{40}$"), x = "The number of chain state", y = "Value")

lags <- seq(0, 100, 1)
autocorr.beta <- cbind.data.frame(autocorr(mc1$chain.beta[,1], lags = lags), 
                                  autocorr(mc1$chain.beta[,nitems/2], lags = lags),
                                  autocorr(mc1$chain.beta[,nitems], lags = lags))
colnames(autocorr.beta) <- c("beta1", paste0("beta", nitems/2), paste0("beta", nitems))
autocorr.beta %>% mutate(index = lags) -> autocorr.beta
autocorr.beta  %>% gather("beta", "value", 1:(dim(autocorr.beta)[2] - 1)) -> autocorr.beta
p12 <- ggplot(data = autocorr.beta, aes(x = index, y = value, color = beta)) + 
  geom_line() + 
  geom_point() +
  geom_abline(aes(slope = 0, intercept = 0)) +
  theme(legend.title = element_blank()) +
  labs(title = TeX("Autocorrelation function for chain $\\beta_{1}, \\beta_{20}, \\beta_{40}$"), y = "Autocorrelation", x = "Lags")
beep()
