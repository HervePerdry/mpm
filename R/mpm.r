mpm <- function(formula, data, rescale.var = TRUE, algorithm = c("scoring", "Nelder-Mead"), control = list()) {
  call <- match.call()
  mf <- model.frame(formula, data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data)
  
  ctr <- list(scoring.maxit = 100, scoring.eps = 1e-12, nm.maxit = 5e4, nm.reltol = 1e-10)
  ctr[ names(control) ] <- control

  X1 <- X[-1,,drop = FALSE]
  if(rescale.var) {
    # center and rescale columns of X to avoid some grossly ill conditionned matrices...
    sc <- apply(X1, 2, sd)
    ce <- colMeans(X1)
    inter <- (sc == 0)
    if(sum(inter) > 1) stop("A problem with the predictors...")
    ce[ inter ] <- 0 # on ne centre pas l'intercept
    sc[ inter ] <- 1 # on ne rescale pas l'intercept 
    # on essaie de garder des coeffs de taille raisonnable dans la matrice de Fisher
    #sc[!inter] <- sc[!inter] * sqrt(nrow(X1)) * sqrt(Y[1]) 
    sc <- sc * sqrt(nrow(X1)) * sqrt(Y[1]) 
  } else {
    sc <- rep(1, ncol(X1))
    ce <- rep(0, ncol(X1)) 
    inter <- rep(FALSE, ncol(X1))
  }
  X1 <- sweep(X1, 2, ce, "-")
  X1 <- sweep(X1, 2, sc, "/")

  # fit model
  algorithm <- match.arg(algorithm)
  
  if(algorithm == "scoring") {
    R <- scoring(Y, X1, max.iter = ctr$scoring.maxit, epsilon = ctr$scoring.eps)
    if(!R$converged) {
      warning("Fisher scoring did not converge, falling back to Nelder-Mead")
      algorithm <- "Nelder-Mead"
    }
  }
  
  if(algorithm == "Nelder-Mead") {
    R <- NM(Y, X1, maxit = ctr$nm.maxit, reltol = ctr$nm.reltol)
    if(!R$converged) 
      warning("Nelder-Mead did not converge")
  }

  # go back to original scale !
  R$beta <- R$beta/sc
  R$inverse.fisher <- R$inverse.fisher / outer(sc, sc)

  # and recompute intercept + its var and covars with others
  R$beta[ inter ] <- R$beta[ inter ] - sum(R$beta*ce)
  z <- -ce; 
  z[ inter ] <- 1
  cv <- R$inv %*% z
  cv[1] <- sum(cv*z)
  R$inverse.fisher[,inter] <- cv
  R$inverse.fisher[inter,] <- cv

  # compute sds, z score, p value
  sds <- sqrt(diag(R$inverse.fisher) * R$phi)
  z <- R$beta / sds
  p.val <- pchisq(z**2, df = 1, lower.tail = FALSE)
  R$coeff <- cbind(beta = R$beta, sd = sds, z.score = z, p.val = p.val)
  colnames(R$coeff) <- c("Estimate", "Std error", "z value", "Pr(>|z|)")
  R$beta <- NULL
  R$null.deviance <- mpm.null.deviance(Y)
  R$algorithm <- algorithm
  R$call <- call
  class(R) <- "mpm"
  R
}
