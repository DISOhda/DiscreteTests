#' @title
#' Poisson Test
#'
#' @description
#' `poisson.test.pv` performs an exact binomial test or a normal
#' approximation about the rate parameter of a Poisson distribution. It is a
#' vectorized version of `poisson.test` that calculates only p-values.
#' Multiple testing scenarios can be passed at once. For two-sided hypotheses,
#' various exact p-value methods are available.
#'
#' @param x              integer vector giving the number of events.
#' @param lambda         hypothesized rate(s).
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x` and `lambda` are vectorised. They are replicated
#' automatically to have the same lengths. This allows multiple hypotheses under
#' the same conditions to be specified simultaneously.
#'
#' Since the Poisson distribution itself has an infinite support, so have the
#' p-values of exact Poisson tests. Thus supports only contain p-values that are
#' not rounded off to 0.
#'
#' @template details_two_sided
#'
#' @template return
#'
#' @seealso
#' [binom.test.pv()], [stats::poisson.test()]
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
#' lambda <- c(3, 2, 1)
#'
#' # Construction of exact two-sided p-values ("blaker") and their supports
#' results.ex  <- poisson.test.pv(k, lambda, ts.method = "blaker")
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_scenario_supports()
#'
#' # Construction of approximate one-sided p-values ("less") and their supports
#' results.ap  <- poisson.test.pv(k, lambda, "less", exact = FALSE)
#' raw.pvalues <- results.ap$get_pvalues()
#' pCDFlist    <- results.ap$get_scenario_supports()
#'
#' @importFrom stats pnorm qnorm
#' @importFrom checkmate assert_atomic_vector assert_integerish assert_numeric
#' @export
poisson.test.pv <- function(
  x,
  lambda = 1,
  alternative = c("two.sided", "less", "greater"),
  ts.method = c("minlike", "blaker", "absdist", "central"),
  exact = TRUE,
  correct = TRUE,
  simple.output = FALSE
) {
  # plausibility checks of input parameters
  len.x <- length(x)
  len.l <- length(lambda)
  len.g <- max(len.x, len.l)

  assert_atomic_vector(x, any.missing = FALSE, min.len = 1)
  assert_integerish(x, lower = 0)
  x <- round(x)
  if(len.x < len.g)
    x <- rep_len(x, len.g)

  assert_atomic_vector(lambda, any.missing = FALSE, min.len = 1)
  assert_numeric(lambda, lower = 0, finite = TRUE)
  if(len.l < len.g)
    lambda <- rep_len(lambda, len.g)

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if(exact && alternative == "two.sided")
    alternative <- match.arg(ts.method,
      c("minlike", "blaker", "absdist", "central")
    )

  lambda.u <- unique(lambda)
  len.u <- length(lambda.u)

  res <- numeric(len.g)
  if(!simple.output) {
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  for(i in 1:len.u) {
    idx <- which(lambda.u[i] == lambda)
    N <- max(x[idx])

    if(exact) {
      # generate all probabilities under current lambda
      d <- generate.poisson.probs(lambda.u[i])
      # difference between number of probabilities and largest desired x-value
      len.diff <- N - length(d) + 1
      # add 0-probabilities, if necessary
      if(len.diff > 0) d <- c(d, rep(0, len.diff))
      # compute p-value support
      pv.supp <- pmin(1,
        switch(alternative,
          less    = cumsum(d),
          greater = rev(cumsum(rev(d))),
          minlike = ts.pv(d, d),
          blaker  = ts.pv(pmin(cumsum(d), rev(cumsum(rev(d)))), d),
          absdist = ts.pv(abs(seq_along(d) - 1 - lambda.u[i]), d,
            decreasing = TRUE
          ),
          central = 2 * pmin(cumsum(d), rev(cumsum(rev(d))))
        )
      )
    } else {
      # find quantile with (approx.) smallest probability > 0
      q <- -floor(qnorm(2^-1074, -lambda.u[i], sqrt(lambda.u[i])))
      # greatest observation
      M <- max(q, N)
      # compute p-value support
      if(lambda.u[i] == 0){
        pv.supp <- switch(alternative,
          less      = rep(1, M + 1),
          greater   = c(1, rep(0, M)),
          two.sided = c(1, rep(0, M))
        )
      } else {
        # possible observations (minus expectation)
        z <- 0:M - lambda.u[i]
        # standard deviation
        std <- sqrt(lambda.u[i])
        # compute p-value support
        pv.supp <- switch(alternative,
          less      = rev(c(1, pnorm(rev(z)[-1] + correct * 0.5, 0, std))),
          greater   = c(1, pnorm(z[-1] - correct * 0.5, 0, std, FALSE)),
          two.sided = pmin(1, 2 * pnorm(-abs(z) + correct * 0.5, 0, std))
        )
      }
    }

    res[idx] <- pv.supp[x[idx] + 1]
    if(!simple.output) {
      supports[[i]] <- sort(unique(pv.supp[pv.supp > 0]))
      indices[[i]]  <- idx
    }
  }

  out <- if(!simple.output) {
    DiscreteTestResults$new(
      test_name = ifelse(exact, "Exact Poisson test",
        paste0("Normal-approximated Poisson test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(`number of events` = x,
        parameters = list(`event rate` = lambda.u)
      ),
      alternative = alternative,
      p_values = res,
      scenario_supports = supports,
      scenario_indices = indices,
      data_name = sapply(match.call(), deparse1)["x"]
    )
  } else res

  return(out)
}
