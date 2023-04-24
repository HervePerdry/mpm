# Fisher scoring 
scoring <- function(Y, X, beta0, max.iter, epsilon) {
  n <- nrow(X)
  p <- ncol(X)
  if(length(Y) != (n+1)) 
    stop("Y doit avoir une observation de plus qu'il n'y a de covariables")

  Y1 <- Y[-(n+1)]
  Y2 <- Y[-1]

  rho.sat <- Y2/Y1
  L0 <- sum(Y2*log(rho.sat) - rho.sat*Y1)

  if(missing(beta0)) { 
    # 'clever' starting point, assuming 1st col is intercept
    rho0 <- sum(Y2)/sum(Y1)
    beta0 <- c(log(rho0), numeric(p-1))/mean(X[,1])
  }
  beta <- beta0

  rho <- exp(as.vector(X %*% beta)); 
  Dev.old <- 2*L0 - 2*sum(Y2*log(rho) - rho*Y1)

  cv <- FALSE # convergence

  for(i in 1:max.iter) {
    U <- score(beta, Y, X, rho)
    I <- fisher.information(beta, Y, X, rho)
    beta <- beta + solve(I, U)

    rho <- exp(as.vector(X %*% beta)); 
    Dev <- 2*L0 - 2*sum(Y2*log(rho) - rho*Y1) 
    if( abs(Dev - Dev.old)/(abs(Dev) + 0.1) < epsilon ) {
      cv <- TRUE
      break
    }
    if( Dev > Dev.old ) {
      cv <- FALSE
      break
    }   
    Dev.old <- Dev
  }

  Y2.hat <- Y1*rho
  residuals <- (Y2 - Y2.hat)/sqrt(Y2.hat)
  phi <- sum(residuals**2) / (n-p)
  list( beta = beta, inverse.fisher = solve(I), phi = phi, pearson.residuals = as.vector(residuals), 
        iterations = i, converged = cv, deviance = Dev)
}
