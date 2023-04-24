mpm.deviance <- function(Y, X, beta, Y1, Y2) {
  if(missing(Y1) | missing(Y2)) {
    n <- nrow(X)
    p <- ncol(X)
    if(length(Y) != (n+1))
      stop("Y doit avoir une observation de plus qu'il n'y a de covariables")

    Y1 <- Y[-(n+1)]
    Y2 <- Y[-1]
  }

  rho.sat <- Y2/Y1
  L0 <- sum(Y2*log(rho.sat) - rho.sat*Y1)

  rho <- exp(as.vector(X %*% beta));
  Dev <- 2*L0 - 2*sum(Y2*log(rho) - rho*Y1)

  # on pourrait utiliser 
  # 2*sum(Y2*log(rho.sat/rho) - (rho.sat - rho)*Y1)
  # est-ce plus stable ?

  Dev
}
