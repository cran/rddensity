################################################################################
#' @title Manipulation Testing Using Local Polynomial Density Estimation
#'
#' @description \code{rddensity} implements manipulation testing procedures using the local
#'   polynomial density estimator proposed in Cattaneo, Jansson and Ma (2019).
#'   For a review on manipulation testing see McCrary (2008).
#'
#' Companion command: \code{\link{rdbwdensity}} for data-driven bandwidth selection, and
#'   \code{\link{rdplotdensity}} for density plot.
#'   A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).
#'
#' Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described at \url{https://sites.google.com/site/rdpackages}.
#'
#' @param X Numeric vector or one dimensional matrix / data frame, the running variable.
#' @param c Numeric, specifies the threshold or cutoff value in the support of \code{X},
#'   which determes the two samples (e.g., control and treatment units in RD settings).  Default
#'   is \code{0}.
#' @param p Integer, specifies the order of the local-polynomial used to construct the density
#'   point estimators.  Default is \code{2} (local quadratic approximation).
#' @param q Integer, specifies the order of the local-polynomial used to construct the
#'   bias-corrected density point estimators.  Default is \code{p+1} (local cubic approximation).
#' @param kernel String, specifies the kernel function used to construct the local-polynomial
#'   estimator(s). Options are: \code{"triangular"}, \code{"epanechnikov"}, and \code{"uniform"}. Default is
#'   \code{"triangular"}.
#' @param fitselect String, specifies whether restrictions should be imposed. Options are:
#'   \code{"unrestricted"} for density estimation without any restrictions (two-sample, unrestricted
#'   inference). This is the default option. \code{"restricted"} for density estimation assuming
#'   equal c.d.f. and higher-order derivatives.
#' @param h Numeric, specifies the bandwidth used to construct the density estimators on the
#'   two sides of the cutoff. If not specified, the bandwidth is computed by the companion
#'   command \code{\link{rdbwdensity}}. If two bandwidths are specified, the first bandwidth
#'   is used for the data below the cutoff and the second bandwidth is used for the data above the cutoff.
#' @param bwselect String, specifies the bandwidth selection procedure to be used. Options are:
#'   \code{"each"} bandwidth selection based on MSE of each density separately (two distinct bandwidths).
#'   \code{"diff"} bandwidth selection based on MSE of difference of densities (one common bandwidth).
#'   \code{"sum"} bandwidth selection based on MSE of sum of densities (one common bandwidth).
#'   \code{"comb"} (this is the default option) bandwidth is selected as a combination of the alternatives above: for \code{fitselect = "unrestricted"},
#'   it selects \code{median(each,diff,sum)}; for \code{fitselect = "restricted"}, it selects \code{min(diff,sum)}.
#' @param vce String, specifies the procedure used to compute the variance-covariance matrix estimator. Options are:
#'   \code{"plugin"} for asymptotic plug-in standard errors. \code{"jackknife"} for jackknife standard errors. This
#'   is the default option.
#' @param all Boolean, if specified, rddensity reports two testing procedures (given choices \code{fitselect}
#'    and \code{bwselect}): Conventional test statistic (not valid when using MSE-optimal bandwidth choice).
#'    Robust bias-corrected statistic.
#'
#' @return
#' \item{hat}{\code{left}/\code{right}: density estimate to the left/right of cutoff; \code{diff}: difference in
#'   estimated densities on the two sides of cutoff.}
#' \item{sd_asy}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
#'   cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on
#'   asymptotic formula.)}
#' \item{sd_jk}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
#'   cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on the
#'   jackknife method.)}
#' \item{test}{\code{t_asy}/\code{t_jk}: t-statistic for the density discontinuity test, with standard error
#'   based on asymptotic formula or the jackknife; \code{p_asy}/\code{p_jk}: p-value for the density
#'   discontinuity test, with standard error based on asymptotic formula or the jackknife.}
#' \item{hat_p}{Same as \code{hat}, without bias correction (only available when \code{all=TRUE}).}
#' \item{sd_asy_p}{Same as \code{sd_asy}, without bias correction (only available when \code{all=TRUE}).}
#' \item{sd_jk_p}{Same as \code{sd_jk}, without bias correction (only available when \code{all=TRUE}).}
#' \item{test_p}{Same as \code{test}, without bias correction (only available when \code{all=TRUE}).}
#' \item{N}{\code{full}: full sample size; \code{left}/\code{right}: sample size to the left/right of the cutoff;
#'   \code{eff_left}/\code{eff_right}: effective sample size to the left/right of the cutoff (this depends
#'   on the bandwidth).}
#' \item{h}{\code{left}/\code{right}: bandwidth used to the left/right of the cutoff.}
#' \item{opt}{Collects the options used, includes: \code{fitselect}, \code{kernel}, \code{bwselectl},
#'   \code{bwselect}, \code{hscale}, \code{vce}, \code{c}, \code{p}, \code{q}, \code{all}.
#'   See options for \code{rddensity}.}
#' \item{X_min}{\code{left}/\code{right}: the samllest observation to the left/right of the cutoff.}
#' \item{X_max}{\code{left}/\code{right}: the largest observation to the left/right of the cutoff.}
#'
#' @references
#' M.D. Cattaneo, M. Jansson and X. Ma. (2018). \href{https://sites.google.com/site/rdpackages/rddensity/Cattaneo-Jansson-Ma_2018_Stata.pdf}{Manipulation Testing based on Density Discontinuity}.  \emph{Stata Journal} 18(1): 234-261.
#'
#' M.D. Cattaneo, M. Jansson and X. Ma. (2019).  \href{https://arxiv.org/abs/1811.11512}{Simple Local Polynomial Density Estimators}. \emph{Journal of the American Statistical Association}, forthcoming.
#'
#' J. McCrary. (2008). Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test. \emph{Journal of Econometrics} 142(2): 698-714.
#'
#' @author
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @seealso \code{\link{rdbwdensity}}, \code{\link{rdplotdensity}}
#'
#' @examples
#' # Continuous Density
#' set.seed(42)
#' x <- rnorm(2000, mean = -0.5)
#' summary(rddensity(X = x, vce="jackknife"))
#'
#' # Discontinuous density
#' x[x>0] <- x[x>0] * 2
#' summary(rddensity(X = x, vce="jackknife"))
#'
#' @export
rddensity <- function(X, c=0, p=2, q=0, kernel="", fitselect="", h=c(), bwselect="", vce="", all=FALSE) {

  ################################################################################
  # default values
  ################################################################################
  if (q==0) { q <- p+1 }
  if (kernel == "") { kernel <- "triangular" }
  kernel <- tolower(kernel)
  if (fitselect == "") { fitselect <- "unrestricted" }
  fitselect <- tolower(fitselect)
  if (bwselect == "") { bwselect <- "comb" }
  bwselect <- tolower(bwselect)
  if (vce == "") { vce <- "jackknife" }
  vce <- tolower(vce)
  # end of default values

  ################################################################################
  # bandwidth values
  ################################################################################
  if (length(h) == 0) {
    hl <- hr <- 0
  }
  if (length(h) == 1) {
    hl <- hr <- h
    if (h <= 0) {
      stop("Bandwidth has to be positive.")
    }
  }
  if (length(h) == 2) {
    hl <- h[1]; hr <- h[2]
    if (min(h) <= 0) {
      stop("Bandwidth has to be positive.")
    }
  }
  if (length(h) > 2) {
    stop("No more than two bandwidths are accepted.")
  }

  ################################################################################
  # missing value handling
  ################################################################################
  X <- as.vector(X)
  if (any(is.na(X))) {
    warning(paste(sum(is.na(X)), " missing ", switch((sum(is.na(X))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    X <- X[!is.na(X)]
  }
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
  if (p >= q) { stop("q should be larger than p.") }
  if (kernel!="uniform" & kernel!= "triangular" & kernel!="epanechnikov") { stop("kernel incorrectly specified.") }
  if (fitselect!="unrestricted" & fitselect!="restricted") { stop("fitselect incorrectly specified.") }
  if (fitselect=="restricted" & hl!=hr) { stop("Bandwidths must be equal in restricted model.") }
  if (bwselect!="each" & bwselect!="diff" & bwselect!="sum" & bwselect!="comb") { stop("bwselect incorrectly specified.") }
  if (fitselect=="restricted" & bwselect=="each") { stop("bwselect=each is not available in the restricted model.") }
  if (vce!="plugin" & vce!="jackknife") { stop("vce incorrectly specified.") }
  # end of error handling

  ################################################################################
  # bandwidth selection
  ################################################################################
  if (hl > 0 & hr > 0) { bwselectl <- "mannual"} else { bwselectl <- "estimated"
    out <- rdbwdensity(X=X, c=c, p=p, kernel=kernel, fitselect=fitselect, vce=vce)$h
    if (fitselect=="unrestricted" & bwselect=="each" & hl==0)  hl = out[1,1]
  	if (fitselect=="unrestricted" & bwselect=="each" & hr==0)  hr = out[2,1]
		if (fitselect=="unrestricted" & bwselect=="diff" & hl==0)  hl = out[3,1]
		if (fitselect=="unrestricted" & bwselect=="diff" & hr==0)  hr = out[3,1]
		if (fitselect=="unrestricted" & bwselect=="sum"  & hl==0)  hl = out[4,1]
		if (fitselect=="unrestricted" & bwselect=="sum"  & hr==0)  hr = out[4,1]
		if (fitselect=="unrestricted" & bwselect=="comb" & hl==0)  hl = median(c(out[1,1], out[3,1], out[4,1]))
		if (fitselect=="unrestricted" & bwselect=="comb" & hr==0)  hr = median(c(out[2,1], out[3,1], out[4,1]))

		if (fitselect=="restricted" & bwselect=="diff" & hl==0)  hl = out[3,1]
		if (fitselect=="restricted" & bwselect=="diff" & hr==0)  hr = out[3,1]
		if (fitselect=="restricted" & bwselect=="sum"  & hl==0)  hl = out[4,1]
		if (fitselect=="restricted" & bwselect=="sum"  & hr==0)  hr = out[4,1]
		if (fitselect=="restricted" & bwselect=="comb" & hl==0)  hl = min(c(out[3,1], out[4,1]))
		if (fitselect=="restricted" & bwselect=="comb" & hr==0)  hr = min(c(out[3,1], out[4,1]))
  }
  # end of bandwidth selection

  ################################################################################
  # data trimming
  ################################################################################
  hscale <- 1;
  X <- X - c; Y <- (0:(N-1)) / (N-1);
  Xh <- X[(X>=-1*hl*hscale) & (X<=hr*hscale)]; Yh <- Y[(X>=-1*hl*hscale) & (X<=hr*hscale)]
  Nlh <- sum(Xh < 0); Nrh <- sum(Xh >= 0)
  if (Nlh < 5 | Nrh < 5) { stop("Not enough observations to perform calculation.") }
  Nh <- Nlh + Nrh
  # end of data trimming

  ################################################################################
  # estimation
  ################################################################################
  fV_q <- rddensity_fV(Y=Yh, X=Xh, Nl=Nl, Nr=Nr, Nlh=Nlh, Nrh=Nrh, hl=hl*hscale, hr=hr*hscale, p=q, s=1, kernel=kernel, fitselect=fitselect, vce=vce)
  T_asy <- fV_q[3, 1] / sqrt(fV_q[3, 3]); T_jk <- fV_q[3, 1] / sqrt(fV_q[3, 2])
  p_asy <- pnorm(abs(T_asy), lower.tail=FALSE) * 2; p_jk <- pnorm(abs(T_jk), lower.tail=FALSE) * 2

  if (all) {
    fV_p <- rddensity_fV(Y=Yh, X=Xh, Nl=Nl, Nr=Nr, Nlh=Nlh, Nrh=Nrh, hl=hl*hscale, hr=hr*hscale, p=p, s=1, kernel=kernel, fitselect=fitselect, vce=vce)
    T_asy_p <- fV_p[3, 1] / sqrt(fV_p[3, 3]); T_jk_p <- fV_p[3, 1] / sqrt(fV_p[3, 2])
    p_asy_p <- pnorm(abs(T_asy_p), lower.tail=FALSE) * 2; p_jk_p <- pnorm(abs(T_jk_p), lower.tail=FALSE) * 2

    result <- list( hat=   list(left=fV_q[1,1], right=fV_q[2,1], diff=fV_q[3,1]),
                    sd_asy=list(left=sqrt(fV_q[1,3]), right=sqrt(fV_q[2,3]), diff=sqrt(fV_q[3,3])),
                    sd_jk= list(left=sqrt(fV_q[1,2]), right=sqrt(fV_q[2,2]), diff=sqrt(fV_q[3,2])),
                    test=  list(t_asy=T_asy, t_jk=T_jk, p_asy=p_asy, p_jk=p_jk),

                    hat_p=   list(left=fV_p[1,1], right=fV_p[2,1], diff=fV_p[3,1]),
                    sd_asy_p=list(left=sqrt(fV_p[1,3]), right=sqrt(fV_p[2,3]), diff=sqrt(fV_p[3,3])),
                    sd_jk_p= list(left=sqrt(fV_p[1,2]), right=sqrt(fV_p[2,2]), diff=sqrt(fV_p[3,2])),
                    test_p=  list(t_asy=T_asy_p, t_jk=T_jk_p, p_asy=p_asy_p, p_jk=p_jk_p),

                    N=     list(full=N, left=Nl, right=Nr, eff_left=Nlh, eff_right=Nrh),
                    h=     list(left=hl*hscale, right=hr*hscale),
                    opt=   list(fitselect=fitselect, kernel=kernel, bwselectl=bwselectl, bwselect=bwselect,
                                vce=vce, c=c, p=p, q=q, all=all),
                    X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                    X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c))
  } else {
    result <- list( hat=   list(left=fV_q[1,1], right=fV_q[2,1], diff=fV_q[3,1]),
                    sd_asy=list(left=sqrt(fV_q[1,3]), right=sqrt(fV_q[2,3]), diff=sqrt(fV_q[3,3])),
                    sd_jk= list(left=sqrt(fV_q[1,2]), right=sqrt(fV_q[2,2]), diff=sqrt(fV_q[3,2])),
                    test=  list(t_asy=T_asy, t_jk=T_jk, p_asy=p_asy, p_jk=p_jk),

                    hat_p=   list(left=NA, right=NA, diff=NA),
                    sd_asy_p=list(left=NA, right=NA, diff=NA),
                    sd_jk_p= list(left=NA, right=NA, diff=NA),
                    test_p=  list(t_asy=NA, t_jk=NA, p_asy=NA, p_jk=NA),

                    N=     list(full=N, left=Nl, right=Nr, eff_left=Nlh, eff_right=Nrh),
                    h=     list(left=hl*hscale, right=hr*hscale),
                    opt=   list(fitselect=fitselect, kernel=kernel, bwselectl=bwselectl, bwselect=bwselect,
                                vce=vce, c=c, p=p, q=q, all=all),
                    X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                    X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c))
  }


  class(result) <- "CJMrddensity"
  return(result)
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMrddensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMrddensity <- function(object, ...) {
  x <- object
  cat("\nRD Manipulation Test using local polynomial density estimation.\n")
  cat("\n")

  cat(paste(format("Number of obs =", width=22), toString(x$N$full), sep="")); cat("\n")
  cat(paste(format("Model =", width=22), x$opt$fitselect, sep="")); cat("\n")
  cat(paste(format("Kernel =", width=22), x$opt$kernel, sep="")); cat("\n")
  if (x$opt$bwselectl!="mannual") {
    cat(paste(format("BW method =", width=22), x$opt$bwselect, sep="")); cat("\n")
  } else {
    cat(paste(format("BW method =", width=22), x$opt$bwselectl, sep="")); cat("\n")
  }
  cat(paste(format("VCE method =", width=22), x$opt$vce, sep="")); cat("\n")
  cat("\n")

  cat(paste(format(paste("Cutoff c = ", toString(round(x$opt$c, 3)), sep=""), width=22), format("Left of c", width=20), format("Right of c", width=20), sep="")); cat("\n")
  cat(paste(format("Number of obs", width=22), format(toString(x$N$left), width=20), format(toString(x$N$right), width=20), sep="")); cat("\n")
  cat(paste(format("Eff. Number of obs", width=22), format(toString(x$N$eff_left), width=20), format(toString(x$N$eff_right), width=20), sep="")); cat("\n")
  #cat(paste(format("Min Running var.", width=22), format(toString(round(x$X_min$left, 3)), width=20), format(toString(round(x$X_min$right, 3)), width=20), sep="")); cat("\n")
  #cat(paste(format("Max Running var.", width=22), format(toString(round(x$X_max$left, 3)), width=20), format(toString(round(x$X_max$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Order est. (p)", width=22), format(toString(x$opt$p), width=20), format(toString(x$opt$p), width=20), sep="")); cat("\n")
  cat(paste(format("Order bias (q)", width=22), format(toString(x$opt$q), width=20), format(toString(x$opt$q), width=20), sep="")); cat("\n")
  cat(paste(format("BW est. (h)", width=22), format(toString(round(x$h$left, 3)), width=20), format(toString(round(x$h$right, 3)), width=20), sep="")); cat("\n")
  cat("\n")

  cat(paste(format("Method", width=22), format("T", width=20), format("P > |T|", width=20), sep="")); cat("\n")
  if (x$opt$all) {
    if (x$opt$vce == "plugin") {
      cat(paste(format("Conventional", width=22), format(toString(round(x$test_p$t_asy, 4)), width=20), format(toString(round(x$test_p$p_asy, 4)), width=20), sep="")); cat("\n")
    } else {
      cat(paste(format("Conventional", width=22), format(toString(round(x$test_p$t_jk, 4)), width=20), format(toString(round(x$test_p$p_jk, 4)), width=20), sep="")); cat("\n")
    }
  }
  if (x$opt$vce == "plugin") {
    cat(paste(format("Robust", width=22), format(toString(round(x$test$t_asy, 4)), width=20), format(toString(round(x$test$p_asy, 4)), width=20), sep="")); cat("\n")
  } else {
    cat(paste(format("Robust", width=22), format(toString(round(x$test$t_jk, 4)), width=20), format(toString(round(x$test$p_jk, 4)), width=20), sep="")); cat("\n")
  }
  cat("\n")
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMrddensity} objects.
#'
#' @keywords internal
#' @export
print.CJMrddensity <- function(x, ...) {
  cat("Call:\n")
  cat("rddensity.\n")
  cat("Sample size:\ ", x$N$full, ". ", "Cutoff: ", x$opt$c, ".\n", sep="")
  cat("Model:\ ", x$opt$fitselect, ". ", "Kernel: ", x$opt$kernel, ". ", "VCE: ", x$opt$vce, sep="")
}
