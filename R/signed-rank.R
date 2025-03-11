#' @name signed_rank_test_pv
#'
#' @title
#' Wilcoxon two-sample signed-rank test
#'
#' @description
#' `signed_rank_test_pv()` performs an exact or approximate Wilcoxon signed-rank
#' test about the differences between two paired groups when the data is not
#' necessarily normally distributed. In contrast to [`stats::wilcox.test()`], it
#' is vectorised and only calculates *p*-values. Furthermore, it is capable of
#' returning the discrete *p*-value supports, i.e. all observable *p*-values
#' under a null hypothesis. Multiple tests can be evaluated simultaneously.
#'
#' @param x,y           numerical vectors forming the samples to be tested or
#'                      lists of numerical vectors for multiple samples.
#' @param d             numerical vector of hypothesised differences(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' The parameters `x`, `y`, `d` and `alternative` are vectorised. If `x` and `y`
#' are lists, they are replicated automatically to have the same lengths. In
#' case `x` or `y` are not lists, they are added to new ones, which are then
#' replicated to the appropriate lengths. This allows multiple hypotheses to be
#' tested simultaneously.
#'
#' In the presence of ties or differences of `x` and `y` that are equal to `d`,
#' computation of exact *p*-values is not possible. This also applies, if the
#' sample size is greater than 1,038, because [`stats::dsignrank`] then produces
#' `Inf`s. Therefore, `exact` is ignored in these cases and *p*-values of the
#' respective test settings are calculated by a normal approximation.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`], [`wilcox_single_test_pv()`],
#' [`mann_whitney_test_pv()`]
#'
#' @references
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. pp. 40-55. \doi{10.1002/9781119196037}
#'
#' @examples
#' # Constructing
#' set.seed(1)
#' r1 <- rnorm(100)
#' r2 <- rnorm(100, 1)
#'
#' # Computation of exact two-sided p-values and their supports
#' results_ex  <- mann_whitney_test_pv(r1, r2)
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- mann_whitney_test_pv(r1, r2, alternative = "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom stats pwilcox
#' @importFrom checkmate qassert qassertr
#' @export
signed_rank_test_pv <- function(
  x,
  y,
  d = 0,
  alternative = "two.sided",
  exact = TRUE,
  correct = TRUE,
  digits_rank = Inf,
  simple_output = FALSE
) {
  # plausibility checks of input parameters
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  qassert(y, c("N+", "L+"))
  if(!is.list(y)) y <- list(y) else qassertr(y, "N+")
  len_y <- length(y)

  qassert(d, "N+()")
  len_d <- length(d)

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
  }

  # replicate inputs to same length
  len_g <- max(len_x, len_y, len_d, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_y < len_g) y <- rep_len(y, len_g)
  if(len_d < len_g) d <- rep_len(d, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  for(i in seq_len(len_g))
    x[[i]] <- x[[i]] - y[[i]]

  res <- wilcox_single_test_pv(
    x, d, alternative, exact, correct, digits_rank, simple_output
  )

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon's signed-rank test",
      inputs = list(
        observations = list(x, y),
        nullvalues = data.frame(`location shift` = d, check.names = FALSE),
        parameters = res$get_inputs()$parameters
      ),
      statistics = res$get_statistics(),
      p_values = res$get_pvalues(),
      pvalue_supports = res$get_pvalue_supports(unique = TRUE),
      support_indices = res$get_support_indices(),
      data_name = paste(dnames["x"], "and", dnames["y"])
    )
  }

  return(out)
}
