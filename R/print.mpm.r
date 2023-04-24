print.mpm <- function(x) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Dispersion parameter phi =", x$phi, "\n")
    cat("Deviance =", sprintf("%.2f", x$deviance))
    cat(" (Null deviance =", sprintf("%.2f", x$null.deviance), ")\n\n")
    cat("Coefficients:\n")
    printCoefmat(x$coeff, P.values  = TRUE, has.Pvalue = TRUE)
    cat("\n")
    if(x$converged) 
      cat("Model fitted with", x$algorithm, "algorithm,", x$iterations, "iterations\n")
    else
      cat("*** MODEL FITTING DID NOT CONVERGE ***\n")
    invisible(x)
}
