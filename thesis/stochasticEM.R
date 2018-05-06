EM <- function(betas.s, thetas.s, discrete.space, npersons, nitems, ndiscrete, Y) { 
  old.likelihood <- 1
  new.likelihood <- 0
  iter_count <- 0
  # E STEP  
  for (iter_count in 1:10000) {
    persons.s <- sampleDist(npersons, thetas.s)
    Y.s <- t(rasch.modelling(persons.s, betas.s))
    n <- sapply(discrete.space,  function(x) sum(persons.s == x))
    r <- sapply(discrete.space, function(x)  {
      if (length(which(persons.s == x)) == 1) {
        Y.s[, which(persons.s == x)]
      } else {
        rowSums(Y.s[, which(persons.s == x)])
      }
    })
    # M STEP
    thetas.s <- n/npersons
    opt <- optim(betas.s, 
                 log.likelihood, 
                 NULL, 
                 n, r,
                 method = "Nelder-Mead", 
                 control = c(maxit = 100, reltol = 1e-8, alpha = 1.0, 
                             beta = 0.5, gamma = 2.0, fnscale = -1)
    )
    betas.s <- opt$par
    old.likelihood <- new.likelihood
    new.likelihood <- opt$value
  }
  print.summary(iter_count, new.likelihood, items, betas.s, thetas.s)
  return(list("beta" = betas.s, "theta" = thetas.s))
}

persons <- sampleDist(npersons)
items <- rnorm(nitems, mean = 3, sd = 1)
rmdata <- rasch.modelling(persons, items)
Y <- t(rmdata)
thetas.s <- rep(1/ndiscrete, ndiscrete) 
betas.s <- rep(3, nitems)

result <- EM(betas.s, thetas.s, discrete.space, npersons, nitems, ndiscrete, Y)
