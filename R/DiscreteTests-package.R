#' @title
#' Vectorised Computation of P-Values and Their Supports for Several Discrete
#' Statistical Tests
#'
#' @docType package
#' @import Rcpp
#' @useDynLib DiscreteTests
#' @name DiscreteTests-package
#'
#' @description
#' This package provides vectorised functions for computing p-values of various
#' discrete statistical tests. Exact and approximate computation methods are
#' provided. For exact p-values, several procedures of determining two-sided
#' p-values are included.
#'
#' Additionally, these functions are capable of returning the discrete p-value
#' supports, i.e. all observable p-values under a null hypothesis. These
#' supports can be used for multiple testing procedures in the
#' [`DiscreteFDR`][DiscreteFDR::DiscreteFDR-package] and
#' [`FDX`][FDX::FDX-package] packages.
#'
#' @references
#' Agresti, A. (2002). *Categorical data analysis* (2nd ed.). New York: John
#'   Wiley & Sons. \doi{10.1002/0471249688}
#'
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. *Canadian Journal of Statistics*,
#'   **28**(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Dixon, W. J. and Mood, A. M. (1946). The statistical sign test.
#'   *Journal of the American Statistical Association*, **41**(236),
#'   pp. 557–566. \doi{10.1080/01621459.1946.10501898}
#'
#' Fisher, R. A. (1935). The logic of inductive inference.
#'   *Journal of the Royal Statistical Society Series A*, **98**, pp.
#'   39–54. \doi{10.2307/2342435}
#'
#' Good, P. (2000). *Permutation Tests: A Practical Guide to Resampling
#'   Methods for Testing Hypotheses*. Second Edition. New York: Springer.
#'   \doi{10.1007/978-1-4757-3235-1}
#'
#' Hirji, K. F. (2006). *Exact analysis of discrete data*. New York: Chapman
#'   and Hall/CRC. pp. 55-83. \doi{10.1201/9781420036190}
#'
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. \doi{10.1002/9781119196037}
#'
#' Mann, H. D. & Whitney, D. R. (1947). On a Test of Whether one of Two Random
#'   Variables is Stochastically Larger than the Other. *Ann. Math. Statist.*,
#'   *18*(1), pp. 50-60. \doi{10.1214/aoms/1177730491}
#'
#' Pitman, E. J. G. (1937). Significance tests which may be applied to
#'   samples from any populations. *Supplement to the Journal of the Royal
#'   Statistical Society*, **4**(1), pp. 119–130. \doi{10.2307/2984124}
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
