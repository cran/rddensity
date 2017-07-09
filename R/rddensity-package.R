################################################################################
#' @title rddensity: Manipulation Testing Using Local-Polynomial Density Estimation
#'
#' @description Density discontinuity test (a.k.a. manipulation test) is commonly employed in
#'   regression discontinuity designs and other treatment effect settings to detect whether there is evidence suggesting
#'   perfect self-selection (manipulation) around a cutoff where a treatment/policy
#'   assignment changes.
#'
#' This package provides tools for conducting the aforementioned statistical
#'   test: \code{\link{rddensity}} to construct local polynomial based density
#'   discontinuity test given a prespecified cutoff, \code{\link{rdbwdensity}} to
#'   perform data-driven bandwidth selection, and \code{\link{rdplotdensity}} to construct density plot near the cutoff.
#'   For a review on manipulation testing see McCrary (2008).
#'
#' For more details, and related \code{Stata} and \code{R} packages
#'   useful for analysis of RD designs, visit \url{https://sites.google.com/site/rdpackages}.
#'
#' @author
#' Matias D. Cattaneo, University of Michigan.  \email{cattaneo@umich.edu}.
#'
#' Michael Jansson, University of California, Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of Michigan. \email{xinweima@umich.edu}.
#'
#' @references
#' M.D. Cattaneo, B. Frandsen, and R. Titiunik. (2015).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate}. \emph{Journal of Causal Inference} 3(1): 1-24.
#'
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017a).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Simple Local Polynomial Density Estimators}. Working Paper, University of Michigan.
#'
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017b). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_Stata.pdf}{rddensity: Manipulation Testing based on Density Discontinuity}. Working Paper, University of Michigan.
#'
#' J. McCrary. (2008). Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test. \emph{Journal of Econometrics} 142(2): 698-714.
#'
#' @importFrom graphics text
#' @importFrom stats dnorm
#' @importFrom stats integrate
#' @importFrom stats median
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats quantile
#' @import lpdensity
#' @import ggplot2
#'
#' @aliases rddensity-package
"_PACKAGE"
