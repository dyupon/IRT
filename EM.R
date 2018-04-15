library(MASS)
library(ggplot2)
library(grid)
set.seed(62)

layout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

arrange <- function(..., nrow=NULL, ncol=NULL) {
  dots <- list(...)
  n <- length(dots)
  if (is.null(nrow) & is.null(ncol)) { 
    nrow = floor(n/2); 
    ncol = ceiling(n/nrow)
  }
  if (is.null(nrow)) {
    nrow = ceiling(n/ncol)
  }
  if (is.null(ncol)) {
    ncol = ceiling(n/nrow)
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow,ncol)))
  ii.p <- 1
  for (ii.row in seq(1, nrow)) {
    ii.table.row <- ii.row
    for (ii.col in seq(1, ncol)) {
      ii.table <- ii.p
      if (ii.p > n) 
        break
      print(dots[[ii.table]], vp = layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

#DATA SIMULATION
rasch.modelling <- function(persons, items) {
  n.items <- length(items)
  n.persons <- length(persons)
  fsmat <- outer(persons, items, "-")
  psolve <- exp(fsmat)/(1 + exp(fsmat))
  R <- (matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
  R
}

npersons <- 100
nitems <- 20
ndiscrete <- 2
discrete.space <- c(2, 5) 
prob <- c(0.4, 0.6)
sampleDist = function(n, prob) { 
  sample(x = discrete.space, n, replace = T, prob)
} 

persons <- sampleDist(npersons, prob)
items <- rnorm(nitems, mean = 3, sd = 1)
rmdata <- rasch.modelling(persons, items)
Y <- t(rmdata)


#APPROXIMATIONS
thetas.s <- unique(runif(4*ndiscrete))[1:ndiscrete]
thetas.s <- thetas.s/sum(thetas.s)
betas.s <- rep(3, nitems) 

#UITLS
probability <- function(beta, theta) {
  return(1 / (1 + exp(beta - theta)))
}

log.likelihood <- function(x, n, r) { # equation 16
  P <- outer(x, discrete.space, probability)
  bet <- rep(0, nitems) 
  for (j in (1:nitems)) {
    bet.k <- rep(0, ndiscrete) 
    for (k in (1:ndiscrete)) {
      bet.k[k] <- r[j,k]*log(P)[j,k] + (n[k] - r[j,k])*log(1 - P)[j,k]
    }
    bet[j] <- sum(bet.k)
  }
  res <- sum(bet)
  return(res)
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
    #if (iter_count %% 10 == 0) {
      print("iteration")
      print(iter_count)
      print(betas.s)
      print(thetas.s)
      print(new.likelihood)
    #}
      # p1 <- ggplot(as.data.frame(cumsum(thetas.s)), aes(cumsum(thetas.s))) + stat_ecdf(geom = "step")
      # p2 <- ggplot(as.data.frame(cumsum(prob)), aes(cumsum(prob))) + stat_ecdf(geom = "step")
      # p3 <- ggplot(as.data.frame(betas.s), aes(betas.s)) + geom_density()
      # p4 <- ggplot(as.data.frame(items), aes(items)) + geom_density()
      # arrange(p1,p2,p3,p4, ncol = 2)
  }
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
  c(0.4, 0.6)
  "estimated theta:"
  thetas.s
  return(list("beta" = betas.s, "theta" = thetas.s))
}
