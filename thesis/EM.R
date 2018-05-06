library(MASS)
library(ggplot2)
library(grid)
set.seed(42)
npersons <- 1000
nitems <- 20
ndiscrete <- 2
discrete.space <- c(2, 5) 
true.theta <- c(0.4, 0.6)
sampleDist = function(n) { 
  sample(x = discrete.space, n, replace = T, prob = true.theta)
} 

persons <- sampleDist(npersons)
items <- rnorm(nitems, mean = 3, sd = 1)
rmdata <- rasch.modelling(persons, items)
Y <- t(rmdata)

#APPROXIMATIONS
thetas.s <- rep(1/ndiscrete, ndiscrete) 
betas.s <- rep(3, nitems) 

#UTILS
probability <- function(beta, theta) {
  return(1 / (1 + exp(beta - theta)))
}

log.likelihood <- function(x, n, r, iter_count) {
  P <- outer(x, discrete.space, probability)
  bet <- rep(0, nitems) 
  for (j in (1:nitems)) {
    bet.k <- rep(0, ndiscrete) 
    for (k in (1:ndiscrete)) {
      bet.k[k] <- r[j,k]*log(P)[j,k] + (n[k] - r[j,k])*log(1 - P)[j,k]
    }
    bet[j] <- sum(bet.k)
  }
  sum(bet)
}

print.summary <- function(iter_count, new.likelihood, items, betas.s, thetas.s) {
  "amount of iterations done:"
  iter_count
  "last loglikelihood value"
  new.likelihood
  "true beta:"
  items
  "estimated beta:" 
  betas.s
  "fitdistr beta:"
  fitdistr(betas.s, "normal")
  "true theta:"
  true.theta
  "estimated theta:"
  thetas.s
}
eps <- 10^-6

#P -- matrix w/ probabilities JxK (nitems x ndisrete)
#Y -- response matrix JxN  (nitems x npersons)
#pi_s -- vector of length K (ndiscrete)
#r_s -- matrix JxK (nitems x ndisrete)

EM <- function(betas.s, thetas.s, discrete.space, npersons, nitems, ndiscrete, Y) { 
  old.likelihood <- 1
  new.likelihood <- 0
  iter_count <- 0
  # E STEP  
  while (abs(old.likelihood - new.likelihood) > eps) {
    iter_count <- iter_count + 1
    P <- outer(betas.s, discrete.space, probability)
    prod_j <- exp(colSums(log(P)[, rep(1:ndiscrete, npersons)] * Y[, rep(1:npersons, each = ndiscrete)] + log(1 - P)[, rep(1:ndiscrete, npersons)] * (1 - Y)[, rep(1:npersons, each = ndiscrete)]))
    dim(prod_j) <- c(ndiscrete, npersons)
    prod_j_pi <- prod_j * thetas.s
    prod_j_pi_k <- t(t(prod_j_pi)/rowSums(t(prod_j_pi)))
    n <- rowSums(prod_j_pi_k) #equation 13
    r <- rowSums(prod_j_pi_k[rep(1:ndiscrete, nitems), ] * Y[rep(1:nitems, each = ndiscrete), ])
    r <- matrix(r,nrow = nitems, ncol = ndiscrete, byrow = TRUE) 
    # M STEP
    thetas.s <- n/npersons
    opt <- optim(betas.s, 
                 log.likelihood, 
                 NULL, 
                 n, r, iter_count,
                 method = "Nelder-Mead", 
                 control = c(maxit = 100, reltol = 1e-8, alpha = 1.0, beta = 0.5, gamma = 2.0, fnscale = -1))
    betas.s <- opt$par
    old.likelihood <- new.likelihood
    new.likelihood <- opt$value
  }
  print.summary(iter_count, new.likelihood, items, betas.s, thetas.s)
  return(list("beta" = betas.s, "theta" = thetas.s))
}
result <- EM(betas.s, thetas.s, discrete.space, npersons, nitems, ndiscrete, Y)