#' @details
#' For exact computation, various procedures of determining two-sided p-values
#' are implemented.
#' \describe{
#'   \item{`"minlike"`}{The standard approach in [stats::fisher.test()] and
#'                      [stats::binom.test()]. The probabilities of the
#'                      likelihoods that are equal or less than the observed one
#'                      are summed up. In Hirji (2006), it is referred to as the
#'                      *Probability-based* approach.}
#'   \item{`"blaker"`}{The minima of the observations' lower and upper tail
#'                     probabilities are combined with the opposite tail not
#'                     greater than these minima. More details can be found in
#'                     Blaker (2000) or Hirji (2006), where it is referred to as
#'                     the *Combined Tails* method.}
#'   \item{`"absdist"`}{The probabilities of the absolute distances from the
#'                      expected value that are greater than or equal to the
#'                      observed one are summed up. In Hirji (2006), it is
#'                      referred to as the *Distance from Center* approach.}
#'   \item{`"central"`}{The smaller values of the observations' simply doubles
#'                      the minimum of lower and upper tail probabilities. In
#'                      Hirji (2006), it is referred to as the *Twice the
#'                      Smaller Tail* method.}
#' }
#' For non-exact (i.e. continuous approximation) approaches, `ts.method` is
#' ignored, since all its methods would yield the same p-values. More
#' specifically, they all converge to the doubling approach as in
#' `ts.mthod = "central"`.
