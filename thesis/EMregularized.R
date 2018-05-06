log.likelihood <- function(x, n, r, iter_count) {
  kk <- 0.4
  lambda <- 100
  P <- outer(x, discrete.space, probability)
  bet <- rep(0, nitems) 
  for (j in (1:nitems)) {
    bet.k <- rep(0, ndiscrete) 
    for (k in (1:ndiscrete)) {
      bet.k[k] <- r[j,k]*log(P)[j,k] + (n[k] - r[j,k])*log(1 - P)[j,k]
    }
    bet[j] <- sum(bet.k)
  }
  res <- sum(bet) - kk^iter_count*lambda*sum(x^2)
  return(res)
}

persons <- sampleDist(npersons)
items <- rnorm(nitems, mean = 3, sd = 1)
rmdata <- rasch.modelling(persons, items)
Y <- t(rmdata)
thetas.s <- rep(1/ndiscrete, ndiscrete) 
betas.s <- rep(3, nitems)

result <- EM(betas.s, thetas.s, discrete.space, npersons, nitems, ndiscrete, Y)