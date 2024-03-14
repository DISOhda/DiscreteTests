#' @title
#' Poisson Test
#'
#' @description
#' \code{poisson.test.pv} performs an exact binomial test or a normal
#' approximation about the rate parameter of a Poisson distribution. It is a
#' vectorized version of \code{poisson.test} that calculates only p-values.
#' Multiple testing scenarios can be passed at once. For two-sided hypotheses,
#' various exact p-value methods are available.
#'
#' @param x              integer vector giving the number of events.
#' @param lambda         hypothesized rate(s).
#' @param alternative    indicates the alternative hypothesis and must be one of
#'                       \code{"two.sided"} (the default), \code{"less"} or
#'                       \code{"greater"}.
#' @param ts.method      indicate the two-sided p-value computation method (if
#'                       \code{alternative = two.sided}) and must be one of
#'                       \code{"minlike"} (the default), \code{"blaker"},
#'                       \code{"absdist"} or \code{"central"} (see details).
#'                       Ignored, if \code{exact = FALSE}.
#' @param exact          logical value that indicates whether p-values are to be
#'                       calculated by exact computation (\code{TRUE}; the
#'                       default) or by a continuous approximation.
#' @param correct        logical value that indicates whether Yates' continuity
#'                       correction for continuous approximations is to be
#'                       applied  (\code{TRUE}; the default) or not. Ignored, if
#'                       \code{exact = TRUE}.
#' @param simple.output  logical value that indicates whether the support sets,
#'                       i.e. all attainable p-values of each testing scenario
#'                       are to be returned, too (see below).
#'
#' @details
#' The parameters \code{x} and \code{lambda} are vectorized. They are replicated
#' automatically to have the same lengths. This allows multiple hypotheses under
#' the same conditions to be specified simultaneously.
#'
#' Since the Poisson distribution itself has an infinite support, so have the
#' p-values of exact Poisson tests. Thus supports only contain p-values that are
#' not rounded off to 0.
#'
#' For exact computation, multiple two-sided p-value methods are available.
#' \describe{
#'   \item{\code{"minlike"}}{identical to the standard approach in
#'                           \code{\link[stats]{fisher.test}} and
#'                           \code{\link[stats]{prop.test}}. The probabilities
#'                           of the likelihoods that are equal or less than the
#'                           observed one are summed up. In Hirji (2006), it is
#'                           referred to as the 'Probability-based' approach.}
#'   \item{\code{"blaker"}}{The minima of the observations' lower and upper tail
#'                          probabilities are combined with the opposite tail
#'                          not greater than these minima. It is more closely
#'                          described in Blaker (2000) or Hirji (2006), where it
#'                          is referred to as the 'Combined Tails' method.}
#'   \item{\code{"absdist"}}{The probabilities of the absolute distances from
#'                           the expected value that are greater than or equal
#'                           to the observed one are summed up. In Hirji (2006),
#'                           it is referred to as the 'Distance from Center'
#'                           approach.}
#'   \item{\code{"central"}}{The smaller values of the observations' simply
#'                           doubles the minimum of lower and upper tail
#'                           probabilities. In Hirji (2006), it is referred to
#'                           as the 'Twice the Smaller Tail' method.}
#' }
#' For non-exact (i.e. continuous approximation) approaches, \code{ts.method} is
#' ignored, since all its methods would yield the same p-values. More specific,
#' they are all converge to the doubling approach as in
#' \code{ts.mthod = "central"}.
#'
#' @return
#' If \code{supports = TRUE}, a list is returned:
#' \item{\code{$p.values}}{a vector of computed p-values as described above.}
#' \item{\code{$supports}}{a list of vectors, each containing all attainable
#'                         p-values of the respective scenario.}
#' Otherwise, only the vector of computed p-values is returned.
#'
#' @seealso
#' \code{\link{binom.test.pv}}, \code{\link[stats]{poisson.test}}
#'
#' @references
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. \emph{Canadian Journal of Statistics},
#'   \strong{28}(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Hirji, K. F. (2006). \emph{Exact analysis of discrete data}. New York: Chapman
#'   and Hall/CRC. pp. 55-83
#'
#' @examples
#' # Constructing
#' k <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' lambda <- c(3, 2, 1)
#'
#' # Construction of exact two-sided p-values ("blaker") and their supports
#' results.ex  <- poisson.test.pv(k, lambda, ts.method = "blaker", simple.output = FALSE)
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_support_values()
#'
#' # Construction of approximate one-sided p-values ("less") and their supports
#' results.ap  <- poisson.test.pv(k, lambda, "less", exact = FALSE, simple.output = FALSE)
#' raw.pvalues <- results.ap$get_pvalues()
#' pCDFlist    <- results.ap$get_support_values()
#'
#' @importFrom stats pnorm qnorm
#' @export
poisson.test.pv <- function(x, lambda = 1, alternative = c("two.sided", "less", "greater"), ts.method = c("minlike", "blaker", "absdist", "central"), exact = TRUE, correct = TRUE, simple.output = FALSE){
  len.x <- length(x)
  len.l <- length(lambda)
  len.g <- max(len.x, len.l)

  if(min(len.x, len.l) < 1L)
    stop("Not enough data!")

  if(is.list(x) || !is.vector(x) || !is.numeric(x) || any(is.na(x) | x < 0 | is.infinite(x) | abs(x - (xr <- round(x))) > 1e-14))
    stop("'x' must be a vector of finite and non-negative integers")
  x <- xr
  if(len.x < len.g) x <- rep_len(x, len.g)

  if(is.list(lambda) || !is.vector(lambda) || !is.numeric(lambda) || any(is.na(lambda) | lambda < 0 | is.infinite(lambda)))
    stop("'lambda' must be a vector of finite and non-negative numbers")
  if(len.l < len.g) lambda <- rep_len(lambda, len.g)

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if(exact && alternative == "two.sided")
    alternative <- match.arg(ts.method, c("minlike", "blaker", "absdist", "central"))

  lambda.u <- unique(lambda)
  len.u <- length(lambda.u)

  res <- numeric(len.g)
  if(!simple.output){
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  for(i in 1:len.u){
    idx <- which(lambda.u[i] == lambda)
    N <- max(x[idx])

    if(exact){
      # generate all probabilities under current lambda
      d <- generate.poisson.probs(lambda.u[i])
      # difference between number of probabilities and largest desired x-value
      len.diff <- N - length(d) + 1
      # add 0-probabilities, if necessary
      if(len.diff > 0) d <- c(d, rep(0, len.diff))
      # compute p-value support
      pv.supp <- pmin(1, switch(alternative,
                                less    = cumsum(d),
                                greater = rev(cumsum(rev(d))),
                                minlike = ts.pv(d, d),
                                blaker  = ts.pv(pmin(cumsum(d), rev(cumsum(rev(d)))), d),
                                absdist = ts.pv(abs(seq_along(d) - 1 - lambda.u[i]), d, decreasing = TRUE),
                                central = 2 * pmin(cumsum(d), rev(cumsum(rev(d))))))
    }else{
      # find quantile with (approx.) smallest probability > 0
      q <- -floor(qnorm(2^-1074, -lambda.u[i], sqrt(lambda.u[i])))
      # greatest observation
      M <- max(q, N)
      # compute p-value support
      if(lambda.u[i] == 0){
        pv.supp <- switch(alternative,
                          less      = rep(1, M + 1),
                          greater   = c(1, rep(0, M)),
                          two.sided = c(1, rep(0, M)))
      }else{
        # possible observations (minus expectation)
        z <- 0:M - lambda.u[i]
        # standard deviation
        std <- sqrt(lambda.u[i])
        # compute p-value support
        pv.supp <- switch(alternative,
                          less      = rev(c(1, pnorm(rev(z)[-1] + correct * 0.5, 0, std))),
                          greater   = c(1, pnorm(z[-1] - correct * 0.5, 0, std, FALSE)),
                          two.sided = pmin(1, 2 * pnorm(-abs(z) + correct * 0.5, 0, std)))
      }
    }

    res[idx] <- pv.supp[x[idx] + 1]
    if(!simple.output){
      supports[[i]] <- sort(unique(pv.supp[pv.supp > 0]))
      indices[[i]]  <- idx
    }
  }

  out <- if(!simple.output){
    dname <- sapply(match.call(), deparse1)["x"]
    DiscreteTestResults$new(
      ifelse(exact, "Exact Poisson Test", "Normal-approximated Poisson Test"),
      list(
        `Number of Events` = x,
        `lambda` = lambda.u),
      alternative,
      res,
      supports,
      indices,
      dname
    )
  }else res

  return(out)
}
