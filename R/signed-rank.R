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
#'                      lists of numerical vectors for multiple samples; all
#'                      sample pairs must have the same length.
#' @param shift         numerical vector of hypothesised differences(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' The parameters `x`, `y`, `shift` and `alternative` are vectorised. If `x` and
#' `y` are lists, they are replicated automatically to have the same lengths. In
#' case `x` or `y` are not lists, they are added to new ones, which are then
#' replicated to the appropriate lengths. This allows multiple hypotheses to be
#' tested simultaneously.
#'
#' In the presence of ties or differences of `x` and `y` that are equal to
#' `shift`, computation of exact *p*-values is not possible. Therefore, `exact`
#' is ignored in these cases and *p*-values of the respective test settings are
#' calculated by a normal approximation.
#'
#' By setting `exact = NULL`, exact computation is performed if the differences
#' of both samples in a test setting do not have any ties or zeros and if both
#' sample sizes are lower than or equal to 200.
#'
#' The used test statistics `W` is also known as \eqn{T+} and is defined as the
#' sum of ranks of all strictly positive values of the sample `x`.
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
#' results_ex  <- signed_rank_test_pv(r1, r2)
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- signed_rank_test_pv(r1, r2, alternative = "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom checkmate qassert qassertr
#' @export
signed_rank_test_pv <- function(
  x,
  y,
  shift = 0,
  alternative = "two.sided",
  exact = NULL,
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

  qassert(shift, "N+()")
  len_s <- length(shift)

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
  }

  # replicate inputs to same length
  len_g <- max(len_x, len_y, len_s, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_y < len_g) y <- rep_len(y, len_g)
  if(len_s < len_g) shift <- rep_len(shift, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  for(i in seq_len(len_g)) {
    # check if lengths are equal; stop if they are not
    if(length(x[[i]]) != length(y[[i]]))
      stop('All paired samples must have the same length')
    # compute differences
    x[[i]] <- x[[i]] - y[[i]]
  }

  res <- wilcox_single_test_pv(
    x, shift, alternative, exact, correct, digits_rank, simple_output
  )

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon's signed-rank test",
      inputs = list(
        observations = list(x, y),
        nullvalues = data.frame(`location shift` = shift, check.names = FALSE),
        parameters = res$get_inputs()$parameters
      ),
      statistics = res$get_statistics(),
      p_values = res$get_pvalues(),
      pvalue_supports = res$get_pvalue_supports(unique = TRUE),
      support_indices = res$get_support_indices(),
      data_name = paste(dnames["x"], "and", dnames["y"])
    )
  } else res

  return(out)
}
