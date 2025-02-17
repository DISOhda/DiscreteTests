#' @name binom_test_pv
#'
#' @title
#' Binomial Tests
#'
#' @description
#' `binom_test_pv()` performs an exact or approximate binomial test about the
#' probability of success in a Bernoulli experiment. In contrast to
#' [`stats::binom.test()`], it is vectorised, only calculates p-values and
#' offers a normal approximation of their computation. Furthermore, it is
#' capable of returning the discrete p-value supports, i.e. all observable
#' p-values under a null hypothesis. Multiple tests can be evaluated
#' simultaneously. In two-sided tests, several procedures of obtaining the
#' respective p-values are implemented.
#'
#' `r lifecycle::badge('deprecated')`\cr
#' **Note**: Please use `binom_test_pv()`! The older `binom.test.pv()` is
#' deprecated in order to migrate to snake case. It will be removed in a future
#' version.
#'
#' @param x   integer vector giving the number of successes.
#' @param n   integer vector giving the number of trials.
#' @param p   numerical vector of hypothesised probabilities of success.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts_method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#'
#' @details
#' The parameters `x`, `n`, `p` and `alternative` are vectorised. They are
#' replicated automatically to have the same lengths. This allows multiple
#' hypotheses to be tested simultaneously.
#'
#' If `p = NULL`, it is tested if the probability of success is 0.5 with
#' the alternative being specified by `alternative`.
#'
#' @template details_two_sided
#'
#' @template return
#'
#' @seealso
#' [`stats::binom.test()`]
#'
#' @references
#' Agresti, A. (2002). *Categorical data analysis* (2nd ed.). New York: John
#'   Wiley & Sons. pp. 14-15. \doi{10.1002/0471249688}
#'
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. *Canadian Journal of Statistics*,
#'   **28**(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Hirji, K. F. (2006). *Exact analysis of discrete data*. New York: Chapman
#'   and Hall/CRC. pp. 55-83. \doi{10.1201/9781420036190}
#'
#' @examples
#' # Constructing
#' k <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' n <- c(18, 12, 10)
#' p <- c(0.5, 0.2, 0.3)
#'
#' # Computation of exact two-sided p-values ("blaker") and their supports
#' results_ex  <- binom_test_pv(k, n, p, ts_method = "blaker")
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- binom_test_pv(k, n, p, "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom stats pnorm
#' @importFrom checkmate assert_integerish qassert
#' @export
binom_test_pv <- function(
  x,
  n,
  p = 0.5,
  alternative = "two.sided",
  ts_method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # catch input values or data names from call
  dnames <- sapply(match.call(), deparse1)
  if(is.na(dnames["p"])) dnames <- c(dnames, p = deparse1(substitute(p)))

  # plausibility checks of input parameters
  qassert(x, "A+")
  x <- assert_integerish(x, lower = 0, min.len = 1, coerce = TRUE)
  len_x <- length(x)

  qassert(n, "A+")
  n <- assert_integerish(n, lower = 0, min.len = 1, coerce = TRUE)
  len_n <- length(n)

  len_p <- length(p)
  if(!len_p) {
    p <- 0.5
  } else {
    qassert(p, "A+")
    qassert(p, "N+[0, 1]")
  }

 if(any(x > n))
    stop("All values of 'x' must not exceed 'n'.")

  qassert(exact, "B1")
  if(!exact) qassert(correct, "B1")

  ts_method <- match.arg(
    ts_method,
    c("minlike", "blaker", "absdist", "central")
  )

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
    if(exact && alternative[i] == "two.sided")
      alternative[i] <- ts_method
  }

  qassert(simple_output, "B1")

  # replicate inputs to same length
  len_g <- max(len_x, len_n, len_p, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_n < len_g) n <- rep_len(n, len_g)
  if(len_p < len_g) p <- rep_len(p, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # determine unique parameter sets
  params <- unique(data.frame(n, p, alternative))
  n_u    <- params$n
  p_u    <- params$p
  alt_u  <- params$alternative
  len_u  <- length(n_u)

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # begin computations
  for(i in seq_len(len_u)) {
    idx_supp <- which(n == n_u[i] & p == p_u[i] & alternative == alt_u[i])

    if(exact) {
      # generate all probabilities under current n and p
      d <- generate_binom_probs(n_u[i], p_u[i])
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alt_u[i],
        probs = generate_binom_probs(n_u[i], p_u[i]),
        expectations = if(alt_u[i] == "absdist")
          abs(0:n_u[i] - n_u[i] * p_u[i])
      )
    } else {
      if(p_u[i] == 0)
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = c(rep(1, n_u[i] + 1)),
          greater   = c(1, rep(0, n_u[i])),
          two.sided = c(1, rep(0, n_u[i]))
        )
      else if(p_u[i] == 1)
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = c(rep(0, n_u[i]), 1),
          greater   = rep(1, n_u[i] + 1),
          two.sided = c(rep(0, n_u[i]), 1)
        )
      else {
        # compute p-value support
        pv_supp <- support_normal(
          alternative = alt_u[i],
          x = 0:n_u[i],
          mean = n_u[i] * p_u[i],
          sd = sqrt(n_u[i] * p_u[i] * (1 - p_u[i])),
          correct = correct
        )
      }
    }

    # store results and support
    res[idx_supp] <- pv_supp[x[idx_supp] + 1]
    if(!simple_output) {
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  out <- if(!simple_output) {
    exact_v <- rep_len(exact, len_g)

    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Exact binomial test",
        paste0(
          "Normal-approximated binomial test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(
        observations = data.frame(
          `number of successes` = x,
          check.names = FALSE
        ),
        nullvalues = data.frame(
          `probability of success` = p,
          check.names = FALSE
        ),
        parameters = data.frame(
          `number of trials` = n,
          alternative = alternative,
          exact = exact_v,
          distribution = ifelse(exact_v, "binomial", "normal"),
          check.names = FALSE
        )
      ),
      statistics = NULL,
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = paste0(dnames["x"], ", ", dnames["n"], " and ", dnames["p"])
    )
  } else res

  return(out)
}

#' @rdname binom_test_pv
#' @export
#' @importFrom lifecycle deprecate_soft
binom.test.pv <- function(
    x,
    n,
    p = 0.5,
    alternative = "two.sided",
    ts.method = "minlike",
    exact = TRUE,
    correct = TRUE,
    simple.output = FALSE
) {
  deprecate_soft("0.2", "binom.test.pv()", "binom_test_pv()")

  binom_test_pv(x, n, p, alternative, ts.method, exact, correct, simple.output)
}
