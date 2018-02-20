library(ggplot2)
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
ndiscrete <- 5 
discrete.space <- c(1,2,3,4,5) 

sampleDist = function(n) { 
  sample(x = discrete.space, n, replace = T, prob = c(0.1, 0.2, 0.3, 0.2, 0.3)) 
} 

persons <- sampleDist(npersons)
items <- rnorm(nitems)
rmdata <- rasch.modeling(persons, items)
Y <- rmdata$data


#APPROXIMATIONS
thetas.s <- rep(1/ndiscrete, ndiscrete) 
betas.s <- rep(0.5, nitems) 

#UITLS
probability <- function(beta, theta) {
  return(1 / (1 + exp(theta - beta)))
}

log.likelihood <- function(x) { # equation 16
  P <- outer(x, discrete.space, probability)
  bet <- rep(0, nitems) 
  for (j in (1:nitems)) {
    bet.k <- rep(0, ndiscrete) 
    for (k in (1:ndiscrete)) {
      bet.k <- r[j,k]*log(P)[j,k] + (n[k] - r[j,k])*log(1 - P)[j,k]
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

old.likelihood <- 1
new.likelihood <- 0
iter_count <- 0
# E STEP  
while (abs(old.likelihood - new.likelihood) > eps) {
  iter_count <- iter_count + 1
  P <- outer(betas.s, discrete.space, probability)
  # сворачиваем по j
  prod_j <- exp(colSums(log(P)[, rep(1:ndiscrete, npersons)] * Y[, rep(1:npersons, each = ndiscrete)] + log(1 - P)[, rep(1:ndiscrete, npersons)] * (1 - Y)[, rep(1:npersons, each = ndiscrete)]))
  # получили матрицу
  dim(prod_j) <- c(ndiscrete, npersons)
  # домножили на пи
  prod_j_pi <- prod_j * thetas.s
  # свернули по k -- нормировали вдоль столбцов на сумму по столбцам
  prod_j_pi_k <- t(t(prod_j_pi)/rowSums(t(prod_j_pi)))
  # свернули по N
  n <- rowSums(prod_j_pi_k) #equation 13
  r <- rowSums(prod_j_pi_k[rep(1:ndiscrete, nitems), ] * Y[rep(1:nitems, each = ndiscrete), ])
  r <- matrix(r,nrow = nitems, ncol = ndiscrete, byrow = TRUE) #equation 14
  
  # M STEP
  thetas.s <- n/npersons
  opt <- optim(betas.s, 
               log.likelihood, 
               NULL, 
               method = "Nelder-Mead", 
               control = c(maxit = 100, reltol = 1e-8, alpha = 1.0, beta = 0.5, gamma = 2.0, fnscale = -1))
  betas.s <- opt$par
  old.likelihood <- new.likelihood
  new.likelihood <- opt$value
  if (iter_count %% 10 == 0) {
    print("iteration")
    print(iter_count)
    print(betas.s)
    print(thetas.s)
    print(new.likelihood)
  }
}
ggplot(as.data.frame(betas.s), aes(betas.s)) + stat_ecdf(geom = "step")
ggplot(as.data.frame(rmdata$items), aes(rmdata$items)) + stat_ecdf(geom = "step")
ggplot(as.data.frame(thetas.s), aes(thetas.s)) + stat_ecdf(geom = "step")
ggplot(as.data.frame(rmdata$persons), aes(rmdata$persons)) + stat_ecdf(geom = "step")
ggplot(as.data.frame(betas.s), aes(betas.s)) + geom_density()
ggplot(as.data.frame(rmdata$items), aes(rmdata$items)) + geom_density()
