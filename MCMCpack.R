library(MCMCpack)
set.seed(62)
#DATA SIMULATION
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
items <- rnorm(nitems, 0, 2)
persons <- rnorm(npersons)
rmdata <- sim.rasch(persons, items)

mc1 <- MCMCirt1d(rmdata, beta.start = 1, b0.beta = 1, B0.beta = 0, ab0 = c(0, 0), AB0 = c(1, 0), store.item = TRUE)
summary(mc1)
#autocorr(mc1)

mc2 <- MCMCirt1d(rmdata, beta.start = 1, b0.beta = 1, B0.beta = 0, store.item = TRUE)
summary(mc2)
#autocorr(mc2)

