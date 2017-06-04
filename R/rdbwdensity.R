################################################################################
#' Bandwidth Selection for Manipulation Testing using Local-Polynomial Density Estimation.
#'
#' \code{rdbwdensity} implements several data-driven bandwidth selection methods
#'   useful to construct manipulation testing procedures using the local polynomial
#'   density estimators proposed in Cattaneo, Jansson and Ma (2017a).
#'
#' Companion command: \code{\link{rddensity}} for density discontinuity (manipulation)
#'   testing. A companion \code{Stata} package is described in Cattaneo,
#'   Jansson and Ma (2017b). Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described at \url{https://sites.google.com/site/rdpackages}.
#'
#' @param X Numeric vector or one dimensional matrix / data frame, the running variable.
#' @param c Numeric, specifies the threshold or cutoff value in the support of \code{X},
#'   which determes the two samples (e.g., control and treatment units in RD settings).  Default
#'   is \code{0}.
#' @param p Integer, specifies the order of the local-polynomial used to construct the density
#'   point estimators.  Default is \code{2} (local quadratic approximation).
#' @param kernel String, specifies the kernel function used to construct the local-polynomial
#'   estimator(s). Options are: \code{"triangular"}, \code{"epanechnikov"}, and \code{"uniform"}. Default is
#'   \code{"triangular"}.
#' @param fitselect String, specifies whether restrictions should be imposed. Options are:
#'   \code{"unrestricted"} for density estimation without any restrictions (two-sample, unrestricted
#'   inference). This is the default option. \code{"restricted"} for density estimation assuming
#'   equal c.d.f. and higher-order derivatives.
#' @param vce String, specifies the procedure used to compute the variance-covariance matrix estimator. Options are:
#'   \code{"plugin"} for asymptotic plug-in standard errors. \code{"jackknife"} for jackknife standard errors. This
#'   is the default option.
#'
#' @return
#' \item{h}{Bandwidths for density discontinuity test, left and right to the cutoff, and asymptotic variance and bias.}
#' \item{N}{\code{full}: full sample size; \code{left}/\code{right}: sample size to the left/right of the cutoff.}
#' \item{opt}{Collects the options used, includes: \code{fitselect}, \code{kernel},
#'   \code{vce}, \code{c}, \code{p}. See options for \code{rdbwdensity}.}
#' \item{X_min}{Smallest observations to the left and right of the cutoff.}
#' \item{X_max}{Largest observations to the left and right of the cutoff.}
#'
#' @references
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017a).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Simple Local Regression Distribution Estimators}. Working Paper, University of Michigan.
#'
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017b). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_Stata.pdf}{rddensity: Manipulation Testing based on Density Discontinuity}. Working Paper, University of Michigan.
#'
#' @author
#' Matias D. Cattaneo, University of Michigan.  \email{cattaneo@umich.edu}.
#'
#' Michael Jansson, University of California, Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of Michigan. \email{xinweima@umich.edu}.
#'
#' @seealso \code{\link{rddensity}}
#'
#' @examples
#' set.seed(42); x <- rnorm(2000, mean = -0.5)
#' summary(rdbwdensity(X = x, vce="jackknife"))
#'
#' @export
rdbwdensity <- function(X, c=0, p=2, kernel="", fitselect="", vce="") {

  ################################################################################
  # default values
  ################################################################################
  if (kernel == "") { kernel <- "triangular" }
  kernel <- tolower(kernel)
  if (fitselect == "") { fitselect <- "unrestricted" }
  fitselect <- tolower(fitselect)
  if (vce == "") { vce <- "jackknife" }
  vce <- tolower(vce)
  # end of default values

  ################################################################################
  # sample sizes
  ################################################################################
  X <- sort(X)
  N <- length(X); Nl <- sum(X<c); Nr <- sum(X>=c); Xmin <- min(X); Xmax <- max(X)
  # end of sample sizes

  ################################################################################
  # error handling
  ################################################################################
  if (c <= Xmin | c >= Xmax) { stop("The cutoff should be set within the range of the data.") }
  if (Nl <= 10 | Nr <= 10) { stop("Not enough observations to perform calculations.") }
  if (p!=1 & p!=2 & p!=3 & p!=4 & p!=5 & p!=6 & p!= 7) { stop("p must be an integer between 1 and 7.") }
  if (kernel!="uniform" & kernel!= "triangular" & kernel!="epanechnikov") { stop("kernel incorrectly specified.") }
  if (fitselect!="unrestricted" & fitselect!="restricted") { stop("fitselect incorrectly specified.") }
  if (vce!="plugin" & vce!="jackknife") { stop("vce incorrectly specified.") }
  # end of error handling

  ################################################################################
  # select preliminary bandwidth
  ################################################################################
  X <- X - c; Xmu <- mean(X); Xsd <- sd(X)
  fhatb <- 1 / (rddensity_H(Xmu / Xsd, p+2)^2 * dnorm(Xmu / Xsd))
  fhatc <- 1 / (rddensity_H(Xmu / Xsd, p)^2 * dnorm(Xmu / Xsd))
  # these constants are for uniform kernel
  Cb <- c(25884.444444494150957,3430865.4551236177795,845007948.04262602329,330631733667.03808594,187774809656037.3125,145729502641999264,146013502974449876992)
  Cc <- c(4.8000000000000246914,548.57142857155463389,100800.00000020420703,29558225.458100609481,12896196859.612621307,7890871468221.609375,6467911284037581)
  bn <- ((2*p+1)/4 * fhatb * Cb[p] / N)^(1/(2*p+5))
  cn <- (1/(2*p) * fhatc * Cc[p] / N)^(1/(2*p+1))
  bn <- bn * Xsd; cn <- cn * Xsd
  # end of select preliminary bandwidth

  ################################################################################
  # estimate main bandwidth
  ################################################################################
  Y <- (0:(N-1)) / (N-1)
  Yb <- Y[abs(X) <= bn]; Xb <- X[abs(X) <= bn]; Yc <- Y[abs(X) <= cn]; Xc <- X[abs(X) <= cn]
  Nlb <- sum(Xb < 0); Nrb <- sum(Xb >= 0); Nlc <- sum(Xc < 0); Nrc <- sum(Xc >= 0)

  hn <- matrix(NA, ncol=3, nrow=4)
  colnames(hn) <- c("bw", "variance", "biassq"); rownames(hn) <- c("l", "r", "diff", "sum")
  fV_b <- rddensity_fV(Y=Yb, X=Xb, Nl=Nl, Nr=Nr, Nlh=Nlb, Nrh=Nrb, hl=bn, hr=bn, p=p+2, s=p+1, kernel=kernel, fitselect=fitselect)
  fV_c <- rddensity_fV(Y=Yc, X=Xc, Nl=Nl, Nr=Nr, Nlh=Nlc, Nrh=Nrc, hl=cn, hr=cn, p=p, s=1, kernel=kernel, fitselect=fitselect)

  if (vce == "plugin") { hn[, 2] <- N * cn * fV_c[, 3]  } else { hn[, 2] <- N * cn * fV_c[, 2] }
  if (fitselect == "unrestricted") {
    S <- Sgenerate(p=p, low=0, up=1, kernel=kernel); C <- Cgenerate(k=p+1, p=p, low=0, up=1, kernel=kernel)
    hn[1, 3] <- fV_b[1, 4] * (solve(S) %*% C)[2] * (-1)^p
    hn[2, 3] <- fV_b[2, 4] * (solve(S) %*% C)[2]
    hn[3, 3] <- hn[2, 3] - hn[1, 3]; hn[4, 3] <- hn[2, 3] + hn[1, 3]
  } else {
    Splus <- Splusgenerate(p=p, kernel=kernel); Cplus <- Cplusgenerate(k=p+1, p=p, kernel=kernel); Psi <- Psigenerate(p=p)
    Sinv <- solve(fV_c[2, 1] * Splus + fV_c[1, 1] * Psi%*%Splus%*%Psi)
    C <- fV_b[1, 4] * (fV_c[2, 1] * Cplus + (-1)^(p+1) * fV_c[1, 1] * Psi%*%Cplus)
    temp <- Sinv%*%C
    hn[1, 3] <- temp[2]; hn[2, 3] <- temp[3]; hn[3, 3] <- hn[2, 3] - hn[1, 3]; hn[4, 3] <- hn[2, 3] + hn[1, 3]
  }

  hn[, 3] <- hn[, 3]^2
  hn[, 1] <- (1/(2*p) * hn[, 2] / hn[, 3] / N)^(1/(2*p+1))
  # end of estimate main bandwidth

  for (i in 1:4) {
    if (hn[i, 2] < 0) { hn[i, 1] <- 0; hn[i, 2] <- NA }
    if (is.na(hn[i, 1])) { hn[i, 1] <- 0 }
  }

  result <- list(h=hn, N=list(full=N, left=Nl, right=Nr),
                 opt=list(fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p),
                 X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                 X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c))

  #result <- list(hn=hn, N=N, fitselect=fitselect, kernel=kernel, vce=vce, c=c, Nl=Nl, Nr=Nr, p=p,
  #               X.min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
  #               X.max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c))
  class(result) <- "CJMrdbwdensity"
  return(result)
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMrdbwdensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMrdbwdensity <- function(object, ...) {
  x <- object
  cat("\nBandwidth selection for manipulation testing.\n")
  cat("\n")

  cat(paste(format("Number of obs =", width=20), toString(x$N$full), sep="")); cat("\n")
  cat(paste(format("Model =", width=20), x$opt$fitselect, sep="")); cat("\n")
  cat(paste(format("Kernel =", width=20), x$opt$kernel, sep="")); cat("\n")
  cat(paste(format("VCE method =", width=20), x$opt$vce, sep="")); cat("\n")
  cat("\n")

  cat(paste(format(paste("Cutoff c = ", toString(round(x$opt$c, 3)), sep=""), width=20), format("Left of c", width=20), format("Right of c", width=20), sep="")); cat("\n")
  cat(paste(format("Number of obs", width=20), format(toString(x$N$left), width=20), format(toString(x$N$right), width=20), sep="")); cat("\n")
  cat(paste(format("Min Running var.", width=20), format(toString(round(x$X_min$left, 3)), width=20), format(toString(round(x$X_min$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Max Running var.", width=20), format(toString(round(x$X_max$left, 3)), width=20), format(toString(round(x$X_max$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Order est. (p)", width=20), format(toString(x$opt$p), width=20), format(toString(x$opt$p), width=20), sep="")); cat("\n")
  cat("\n")

  cat(paste(format("Target", width=20), format("Bandwidth", width=20), format("Variance", width=20), format("Bias^2", width=20), sep="")); cat("\n")
  cat(paste(format("left density", width=20), format(toString(round(x$h[1,1], 4)), width=20), format(toString(round(x$h[1,2], 4)), width=20), format(toString(round(x$h[1,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("right density", width=20), format(toString(round(x$h[2,1], 4)), width=20), format(toString(round(x$h[2,2], 4)), width=20), format(toString(round(x$h[2,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("diff. densities", width=20), format(toString(round(x$h[3,1], 4)), width=20), format(toString(round(x$h[3,2], 4)), width=20), format(toString(round(x$h[3,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("sum densities", width=20), format(toString(round(x$h[4,1], 4)), width=20), format(toString(round(x$h[4,2], 4)), width=20), format(toString(round(x$h[4,3], 4)), width=20), sep="")); cat("\n")
  cat("\n")
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMrdbwdensity} objects.
#'
#' @keywords internal
#' @export
print.CJMrdbwdensity <- function(x, ...) {
  cat("Call:\n")
  cat("rdbwdensity\n")
  cat("Sample size:\ ", x$N$full, ". ", "Cutoff: ", x$opt$c, ".\n", sep="")
  cat("Model:\ ", x$opt$fitselect, ". ", "Kernel: ", x$opt$kernel, ". ", "VCE: ", x$opt$vce, sep="")
}







