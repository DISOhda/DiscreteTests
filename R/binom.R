#' @title
#' Binomial Tests
#'
#' @description
#' `binom.test.pv()` performs an exact or approximate binomial test about the
#' probability of success in a Bernoulli experiment. In contrast to
#' [stats::binom.test()], it is vectorised, only calculates p-values and offers
#' a normal approximation of their computation. Furthermore, it is capable of
#' returning the discrete p-value supports, i.e. all observable p-values under a
#' null hypothesis. Multiple tests can be evaluated simultaneously. In two-sided
#' tests, several procedures of obtaining the respective p-values are
#' implemented.
#'
#' @param x   integer vector giving the number of successes.
#' @param n   integer vector giving the number of trials.
#' @param p   numerical vector of hypothesised probabilities of success.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
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
#' [stats::binom.test()]
#'
#' @references
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. *Canadian Journal of Statistics*,
#'   **28**(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Hirji, K. F. (2006). *Exact analysis of discrete data*. New York: Chapman
#'   and Hall/CRC. pp. 55-83
#'
#' @examples
#' # Constructing
#' k <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' n <- c(18, 12, 10)
#' p <- c(0.5, 0.2, 0.3)
#'
#' # Computation of exact two-sided p-values ("blaker") and their supports
#' results.ex  <- binom.test.pv(k, n, p, ts.method = "blaker")
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values (less) and their supports
#' results.ap  <- binom.test.pv(k, n, p, "less", exact = FALSE)
#' raw.pvalues <- results.ap$get_pvalues()
#' pCDFlist    <- results.ap$get_pvalue_supports()
#'
#' @importFrom stats pnorm
#' @importFrom checkmate assert_integerish qassert
#' @export
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
  # plausibility checks of input parameters
  len.x <- length(x)
  len.n <- length(n)
  len.p <- length(p)
  len.a <- length(alternative)
  len.g <- max(len.x, len.n, len.p, len.a)

  qassert(x, "A+")
  x <- assert_integerish(x, lower = 0, min.len = 1, coerce = TRUE)
  if(len.x < len.g) x <- rep_len(x, len.g)

  if(!len.p) {
    p <- rep(0.5, len.g)
  } else {
    qassert(p, "A+")
    qassert(p, "N+[0, 1]")
    if(len.p < len.g) p <- rep_len(p, len.g)
  }

  qassert(n, "A+")
  n <- assert_integerish(n, lower = 0, min.len = 1, coerce = TRUE)
  if(len.n < len.g) n <- rep_len(n, len.g)

  if(any(x > n))
    stop("All values of 'x' must not exceed 'n'.")

  qassert(exact, "B1")
  if(!exact) qassert(correct, "B1")

  ts.method <- match.arg(
    ts.method,
    c("minlike", "blaker", "absdist", "central")
  )

  for(i in seq_len(len.a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
    if(exact && alternative[i] == "two.sided")
      alternative[i] <- ts.method
  }
  if(len.a < len.g) alternative <- rep_len(alternative, len.g)

  qassert(simple.output, "B1")

  # find unique parameter sets
  params <- unique(data.frame(n, p, alternative))
  n.u    <- params$n
  p.u    <- params$p
  alt.u  <- params$alternative
  len.u  <- length(n.u)

  # prepare output
  res <- numeric(len.g)
  if(!simple.output) {
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  # begin computations
  for(i in 1:len.u) {
    idx <- which(n == n.u[i] & p == p.u[i] & alternative == alt.u[i])

    if(exact) {
      # generate all probabilities under current n and p
      d <- generate.binom.probs(n.u[i], p.u[i])
      # compute p-value support
      pv.supp <- support_exact(
        alternative = alt.u[i],
        probs = generate.binom.probs(n.u[i], p.u[i]),
        expectations = if(alt.u[i] == "absdist")
          abs(0:n.u[i] - n.u[i] * p.u[i])
      )
    } else {
      if(p.u[i] == 0)
        pv.supp <- switch(
          EXPR      = alt.u[i],
          less      = c(rep(1, n.u[i] + 1)),
          greater   = c(1, rep(0, n.u[i])),
          two.sided = c(1, rep(0, n.u[i]))
        )
      else if(p.u[i] == 1)
        pv.supp <- switch(
          EXPR      = alt.u[i],
          less      = c(rep(0, n.u[i]), 1),
          greater   = rep(1, n.u[i] + 1),
          two.sided = c(rep(0, n.u[i]), 1)
        )
      else {
        # possible observations (minus expectation)
        z <- 0:n.u[i] - n.u[i] * p.u[i]
        # standard deviation
        std <- sqrt(n.u[i] * p.u[i] * (1 - p.u[i]))
        # compute p-value support
        pv.supp <- switch(
          EXPR      = alt.u[i],
          less      = rev(c(1, pnorm(rev(z)[-1] + correct * 0.5, 0, std))),
          greater   = c(1, pnorm(z[-1] - correct * 0.5, 0, std, FALSE)),
          two.sided = pmin(1, 2 * pnorm(-abs(z) + correct * 0.5, 0, std))
        )
      }
    }

    res[idx] <- pv.supp[x[idx] + 1]
    if(!simple.output) {
      supports[[i]] <- sort(unique(pv.supp))
      indices[[i]]  <- idx
    }
  }

  out <- if(!simple.output) {
    dnames <- sapply(match.call(), deparse1)
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
          check.names = FALSE
        )
      ),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = paste0(dnames["x"], ", ", dnames["n"], " and ", dnames["p"])
    )
  } else res

  return(out)
}
