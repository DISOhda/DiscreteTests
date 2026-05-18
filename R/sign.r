#' @name sign_test_pv
#'
#' @title
#' Sign Tests
#'
#' @description
#' `sign_test_pv()` performs an exact or approximate sign test about the median
#' of a distribution. It supports both one-sample and two-sample paired tests.
#' In contrast to other implementations, it is vectorised, only calculates
#' *p*-values, and offers a normal approximation of their computation.
#' Furthermore, it is capable of returning the discrete *p*-value supports, i.e.
#' all observable *p*-values under a null hypothesis. Multiple tests can be
#' evaluated simultaneously. In two-sided tests, several procedures of obtaining
#' the respective *p*-values are implemented.
#'
#' @param x       numerical vector forming the sample to be tested or a list of
#'                numerical vectors for multiple tests.
#' @param y       numerical vector forming the second sample to be tested or a
#'                list of numerical vectors for multiple tests; if `y = NULL`
#'                (the default), the one-sample version is performed; for
#'                two-sample tests, all sample pairs must have the same length.
#' @param mu      numerical vector of hypothesised median(s) for one-sample
#'                tests or median(s) of differences for two-sample tests.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#'
#' @details
#' For the **one-sample** case, the sign test counts the number of observations
#' in `x` that exceed the hypothesised median `mu`. Observations exactly equal
#' to `mu` (ties) are discarded. Under the null hypothesis, the number of
#' positive signs follows a Binomial(\eqn{n}, 0.5) distribution, where \eqn{n}
#' is the number of non-tied observations.
#'
#' For the **two-sample paired** case (when `y` is supplied), the pairwise
#' differences \eqn{d_i = x_i - y_i - \mu} are computed first. The sign test is
#' then applied to these differences.
#'
#' If `x` is a list, multiple tests are evaluated simultaneously. The parameters
#' `x`, `y` (if supplied and a list), `mu`, and `alternative` are vectorised.
#' They are replicated automatically to have the same lengths. This allows
#' multiple hypotheses to be tested simultaneously.
#'
#' The sign test is a special case of the binomial test. In contrast to
#' [`binom_test_pv()`], `sign_test_pv()` does not allow specifying exact
#' two-sided *p*-value calculation procedures. The reason is that the exact sign
#' test always tests for a probability of 0.5, in which case all these exact
#' two-sided *p*-value computation methods yield exactly the same results.
#'
#' @template return
#'
#' @seealso
#' [binom_test_pv()], [stats::binom.test()]
#'
#' @references
#' Dixon, W. J. and Mood, A. M. (1946). The statistical sign test.
#'   *Journal of the American Statistical Association*, **41**(236),
#'   pp. 557–566. \doi{10.1080/01621459.1946.10501898}
#'
#' @examples
#' # One-sample sign test: test whether median equals 3
#' x <- c(1, 5, 2, 7, 6, 8, 4)
#' results_1s <- sign_test_pv(x, mu = 3)
#' print(results_1s)
#' results_1s$get_pvalues()
#' results_1s$get_pvalue_supports()
#'
#' # Paired two-sample sign test: test whether difference of medians equals 1
#' x2 <- c(6, 8, 2, 4, 5)
#' y2 <- c(3, 5, 4, 2, 6)
#' results_2s <- sign_test_pv(x2, y2, mu = 1)
#' print(results_2s)
#' results_2s$get_pvalues()
#'
#' # Multiple tests simultaneously, one-sided p-values (one-sample, list input)
#' xs <- list(c(1, 5, 2, 7, 6, 8, 4), c(2, 4, 6, 1, 9))
#' results_l <- sign_test_pv(xs, mu = c(3, 5), alternative = "greater")
#' print(results_l)
#' results_l$get_pvalues()
#'
#' # Normal-approximated (one-sample, list input)
#' results_a <- sign_test_pv(xs, mu = c(3, 5), exact = FALSE)
#' print(results_a)
#' results_a$get_pvalues()
#'
#' @importFrom checkmate qassert qassertr
#' @importFrom cli cli_abort cli_warn
#' @export
sign_test_pv <- function(
  x,
  y = NULL,
  mu = 0,
  alternative = "two.sided",
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # plausibility checks of input parameters
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  one_sample <- is.null(y)
  if(!one_sample) {
    qassert(y, c("N+", "L+"))
    if(!is.list(y)) y <- list(y) else qassertr(y, "N+")
  }
  len_y <- length(y)

  qassert(mu, "N+()")
  len_m <- length(mu)

  qassert(exact, c("B1", "0"))
  if(!exact) qassert(correct, "B1")

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      tolower(alternative[i]),
      c("two.sided", "less", "greater")
    )
  }

  qassert(simple_output, "B1")

  # replicate inputs to same length
  len_g <- max(len_x, len_y, len_m, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_m < len_g) mu <- rep_len(mu, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  if(!one_sample) {
    if(len_y < len_g) y <- rep_len(y, len_g)

    for(i in seq_len(len_g)) {
      # check if lengths of all sample pairs are equal; stop if they are not
      if(length(x[[i]]) != length(y[[i]]))
        cli_abort('All paired samples must have the same length')
    }
  }

  # compute sign counts
  s <- integer(len_g)
  n <- integer(len_g)
  for(i in seq_len(len_g)) {
    d    <- if(one_sample) x[[i]] - mu[i] else x[[i]] - y[[i]] - mu[i]
    s[i] <- sum(d > 0)
    n[i] <- sum(d != 0)
  }

  # warn if any test has no non-tied observations
  if(any(n == 0L))
    cli_warn(
      c(
        "Some tests have no non-tied observations; p-values set to 1.",
        i = "Check whether 'x' (and 'y') contain only values equal to 'mu'."
      )
    )

  # determine p-values (and their supports)
  p <- rep(0.5, len_g)
  res <- binom_test_int(s, n, p, alternative, exact, correct, simple_output)

  # create output object
  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = if(one_sample) "One-sample sign test" else
        "Paired two-sample sign test",
      inputs = list(
        observations = if(one_sample) x else list(x, y),
        parameters = data.frame(
          `(non-tied) observations` = n,
          check.names = FALSE
        ),
        nullvalues = if(one_sample) data.frame(median = mu) else
          data.frame(`difference of medians` = mu, check.names = FALSE),
        computation = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            alternative = alternative,
            exact = exact,
            distribution = ifelse(exact, "binomial", "normal"),
            `continuity correction` = if(exact) NA else correct,
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(
        `positive signs` = s,
        #`negative signs` = n - s,
        check.names = FALSE
      ),
      p_values = res$pv,
      pvalue_supports = res$sup,
      support_indices = res$idx,
      data_name = dnames[if(one_sample) "x" else c("x", "y")]
    )
  } else res

  # return results
  return(out)
}
