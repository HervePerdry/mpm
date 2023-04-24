# deviance of the model with only an intercept

mpm.null.deviance <- function(Y, Y1, Y2) {
  if(missing(Y1) | missing(Y2)) {
    n <- length(Y)
    Y1 <- Y[-n]
    Y2 <- Y[-1]
  }

  rho.sat <- Y2/Y1
  L0 <- sum(Y2*log(rho.sat) - rho.sat*Y1)

  # modele avec un rho constant
  rho <- sum(Y2)/sum(Y1)
  Dev <- 2*L0 - 2*sum(Y2*log(rho) - rho*Y1)

  # on pourrait utiliser 
  # 2*sum(Y2*log(rho.sat/rho) - (rho.sat - rho)*Y1)

  Dev
}
