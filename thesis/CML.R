set.seed(42);
library(eRm)
library(ggplot2)
library(MASS)
library(grid)

layout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
arrange <- function(..., nrow=NULL, ncol = NULL) {
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

rasch.modelling <- function(persons, items) {
  n.items <- length(items)
  n.persons <- length(persons)
  fsmat <- outer(persons, items, "-")
  psolve <- exp(fsmat)/(1 + exp(fsmat))
  R <- (matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
  R
}

persons <- rnorm(1000)
items <- rnorm(40, mean = -3, sd = 1)
rmdata <- rasch.modelling(persons, items)
res <- RM(X = rmdata, sum0 = FALSE)
betas <- (res$betapar) 
thetas <- as.vector(person.parameter(res)$thetapar$NAgroup1)
p1 <- ggplot(as.data.frame(thetas), aes(thetas)) + geom_density() +
  ggtitle("Density estimate of thetas")
p2 <- ggplot(as.data.frame(thetas), aes(thetas)) + stat_ecdf(geom = "point") +
  stat_ecdf(geom = "step") + ggtitle("ECDF of thetas")
p3 <- ggplot(as.data.frame(betas), aes(betas)) + geom_density() +
  ggtitle("Density estimate of betas")
p4 <- ggplot(as.data.frame(betas), aes(betas)) + stat_ecdf(geom = "point") +
  stat_ecdf(geom = "step") + labs( title = "ECDF of betas")
cat("Thetas:")
fitdistr(thetas, "normal")
cat("Betas:")  
fitdistr(betas, "normal")
arrange(p1,p2,p3,p4, ncol = 2)
max(thetas)
min(thetas)
max(persons)
min(persons)
