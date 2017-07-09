################################################################################
#' @title Density Plot
#'
#' @description \code{rdplotdensity} constructs density plot near the cutoff. It is based on the
#'   local polynomial density estimators proposed in \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Cattaneo, Jansson and Ma (2017a)}.
#'
#' Companion command: \code{\link{rddensity}} for manipulation testing.
#'   A companion \code{Stata} package is described in \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_Stata.pdf}{Cattaneo, Jansson and Ma (2017b)}.
#'
#' Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described at \url{https://sites.google.com/site/rdpackages}.
#'
#' @param rdd Object returned by \code{\link{rddensity}}
#' @param X Numeric vector or one dimensional matrix / data frame, the running variable.
#' @param plotRange Numeric, specifies the lower and upper bound for density plot. By default it is
#'   three bandwidths around the cutoff.
#' @param plotN Numeric, specifies the number of grid points used for density plot on each side.
#'   If more than one is provided, they will be applied to the two sides accordingly. By default 10
#'   points are used on each side.
#' @param plotGrid String, specifies the position of grid points. Can be either evenly spaced (default, \code{"es"}),
#'   or quantile spaced (\code{"qs"}).
#' @param alpha Numeric scalar between 0 and 1, the significance level for plotting
#'   confidence regions. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param type String, one of \code{"line"} (default), \code{"points"} or \code{"both"}, how
#'   the point estimates are plotted. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param title,xlabel,ylabel Strings, title of the plot and labels for x- and y-axis.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided, they will be applied
#'   to the two sides accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CIshade Numeric, opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CIcol color for confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param legendTitle String, title of legend.
#' @param legendGroups String Vector, group names used in legend.
#'
#' @return
#' \item{Estl, Estr}{Matrices containing estimation results on the two side, with (1) \code{grid} (grid points), (2) \code{bw} (bandwidths), (3) \code{nh}
#'   (effective/local sample sizes), (4) \code{f_p} (point estimates with p-th order local polynomial),
#'   (5) \code{f_q} (point estimates with q-th order local polynomial, only if option \code{q} is nonzero),
#'   (6) \code{se_p} (standard error corresponding to \code{f_p}), and (7) \code{se_q} (standard error
#'   corresponding to \code{f_q}).}
#' \item{Estplot}{A stadnard \code{ggplot} object is returned, hence can be used for further customization.}
#'
#' @references
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017a).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Simple Local Polynomial Density Estimators}. Working Paper, University of Michigan.
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
#' # density being discontinuous
#' set.seed(42)
#' x <- rnorm(2000, mean = -0.5); x[x>0] <- x[x>0] * 2
#' rdd <- rddensity(X = x)
#' plot <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25)
#' @export
rdplotdensity <- function(rdd, X, plotRange = NULL, plotN = 10, plotGrid = c("es", "qs"),
                          alpha = 0.05,
                          type = NULL, CItype = NULL,
                          title = "", xlabel = "", ylabel = "",
                          lty = NULL, lwd = NULL, lcol = NULL,
                          pty = NULL, pwd = NULL, pcol = NULL,
                          CIshade = NULL, CIcol = NULL,
                          legendTitle = NULL, legendGroups = NULL){

  # obtain options from rddensity result
  c       <- rdd$opt$c
  p       <- rdd$opt$p
  q       <- rdd$opt$q
  hl      <- rdd$h$left
  hr      <- rdd$h$right
  kernel  <- rdd$opt$kernel

  # check grid specifications
  if (length(plotRange) == 0) {
    plotRange <- c( max(min(X), c - 3*hl), min(max(X), c + 3 * hr) )
  } else if (length(plotRange) != 2) {
    stop("Plot range incorrectly specified.\n")
  } else if (plotRange[1] >= c | plotRange[2] <= c) {
    stop("Plot range incorrectly specified.\n")
  }

  if (length(plotN) == 0) {
    plotN <- c(10, 10)
  } else if (length(plotN) == 1) {
    plotN <- c(plotN, plotN)
  } else if (length(plotN) > 2) {
    stop("Number of grid points incorrectly specified.\n")
  }
  if (plotN[1] <=1 | plotN[2] <=1) {
    stop("Number of grid points incorrectly specified.\n")
  }

  if (length(plotGrid) == 0) {
    plotGrid <- "es"
  } else {
    plotGrid <- plotGrid[1]
  }
  if (!plotGrid%in%c("es", "qs")) {
    stop("Grid specification invalid.\n")
  }

  # some preparation
  scalel <- (sum(X <= c) - 1) / (length(X) - 1)
  scaler <- (sum(X >= c) - 1) / (length(X) - 1)

  if (plotGrid == "es") {
    gridl <- seq(plotRange[1], c, length.out=plotN[1])
    gridl[plotN[1]] <- c
    gridr <- seq(c, plotRange[2], length.out=plotN[2])
    gridr[1] <- c
  } else {
    gridl <- seq(mean(X <= plotRange[1]), mean(X <= c), length.out=plotN[1])
    gridl <- quantile(X, gridl)
    gridr <- seq(mean(X <= c), mean(X <= plotRange[2]), length.out=plotN[2])
    gridr <- quantile(X, gridr)
    gridl[plotN[1]] <- c
    gridr[1] <- c
  }

  # call lpdensity
  Estl <- lpdensity(data=X[X<=c], grid=gridl, bw=hl, p=p, q=q, v=1, kernel=kernel, scale=scalel)
  Estr <- lpdensity(data=X[X>=c], grid=gridr, bw=hr, p=p, q=q, v=1, kernel=kernel, scale=scaler)

  # call lpdensity.plot
  Estplot <- lpdensity.plot(Estl, Estr, alpha = alpha, type = type, CItype = CItype,
              title = title, xlabel = xlabel, ylabel = ylabel, lty = lty, lwd = lwd,
              lcol = lcol, pty = pty, pwd = pwd, pcol = pcol, CIshade = CIshade,
              CIcol = CIcol, legendTitle = legendTitle, legendGroups = legendGroups) +
    theme(legend.position = "none")

  print(Estplot)

  return(list(Estl=Estl, Estr=Estr, Estplot=Estplot))
}
