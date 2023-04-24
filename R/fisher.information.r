fisher.information <- function(beta, Y, X, rho = exp(X %*% beta)) {
  y0 <- Y[1]
  EW <- y0*cumprod(rho)
  # shorter for t(X) %*% diag(EW) %*% X
  return( crossprod(X, EW * X) )
}
