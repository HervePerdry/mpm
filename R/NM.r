NM <- function(Y, X, beta0) {
  n <- nrow(X)
  p <- ncol(X)
  Y1 <- Y[-(n+1)]
  Y2 <- Y[-1]
  if(missing(beta0)) {
    # 'clever' starting point, assuming 1st col is intercept
    rho0 <- sum(Y2)/sum(Y1)
    beta0 <- c(log(rho0), numeric(p-1))/mean(X[,1])
  }
  beta <- beta0

  o <- optim( beta0, function(b) mpm.deviance(Y, X, b), control = list(reltol = 1e-10))
  beta <- o$par
  I <- fisher.information(beta, Y, X)

  # calcul des résidus
  rho <- exp(as.vector(X %*% beta));
  Y2.hat <- Y1*rho
  residuals <- (Y2 - Y2.hat)/sqrt(Y2.hat)
  phi <- sum(residuals**2) / (n-p)

  list( beta = beta, inverse.fisher = solve(I), phi = phi, pearson.residuals = as.vector(residuals), 
        iterations = as.integer(o$counts[1]), converged = (o$convergence == 0), deviance = o$value)
}
