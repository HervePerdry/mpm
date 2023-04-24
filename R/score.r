score <- function(beta, Y, X, rho = exp(X %*% beta)) {
  n <- nrow(X)
  Y1 <- Y[-(n+1)]
  Y2 <- Y[-1]
  Y2.hat <- Y1*rho
  return( t(X) %*% (Y2 - Y2.hat) )  
}
