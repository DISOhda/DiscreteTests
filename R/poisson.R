#' @title
#' Poisson Test
#'
#' @description
#' `poisson.test.pv()` performs an exact or approximate Poisson test about the
#' rate parameter of a Poisson distribution. In contrast to
#' [stats::poisson.test()], it is vectorised and only calculates p-values.
#' Multiple tests can be evaluated simultaneously. In two-sided tests, several
#' procedures of obtaining the respective p-values are implemented.
#'
#' @param x        integer vector giving the number of events.
#' @param lambda   non-negative numerical vector of hypothesised rate(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x`, `lambda` and `alternative` are vectorised. They are
#' replicated automatically to have the same lengths. This allows multiple
#' hypotheses to be tested simultaneously.
#'
#' Since the Poisson distribution itself has an infinite support, so do the
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
#' # Computation of exact two-sided p-values ("blaker") and their supports
#' results.ex  <- poisson.test.pv(k, lambda, ts.method = "blaker")
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results.ap  <- poisson.test.pv(k, lambda, "less", exact = FALSE)
#' raw.pvalues <- results.ap$get_pvalues()
#' pCDFlist    <- results.ap$get_pvalue_supports()
#'
#' @importFrom stats pnorm qnorm
#' @importFrom checkmate assert_integerish qassert
#' @export
poisson.test.pv <- function(
  x,
  lambda = 1,
  alternative = "two.sided",
  ts.method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple.output = FALSE
) {
  # plausibility checks of input parameters
  len.x <- length(x)
  len.l <- length(lambda)
  len.a <- length(alternative)
  len.g <- max(len.x, len.l, len.a)

  qassert(x, "A+")
  assert_integerish(x, lower = 0)
  x <- round(x)
  if(len.x < len.g) x <- rep_len(x, len.g)

  qassert(x = lambda, "N+[0,)")
  if(len.l < len.g) lambda <- rep_len(lambda, len.g)

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
  params   <- unique(data.frame(lambda, alternative))
  lambda.u <- params$lambda
  alt.u    <- params$alternative
  len.u    <- length(lambda.u)

  # prepare output
  res <- numeric(len.g)
  if(!simple.output) {
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  # begin computations
  for(i in 1:len.u) {
    idx <- which(lambda.u[i] == lambda & alt.u[i] == alternative)
    N <- max(x[idx])

    if(exact) {
      # generate all probabilities under current lambda
      d <- generate.poisson.probs(lambda.u[i])
      # difference between number of probabilities and largest desired x-value
      len.diff <- N - length(d) + 1
      # add 0-probabilities, if necessary
      if(len.diff > 0) d <- c(d, rep(0, len.diff))
      # compute p-value support
      pv.supp <- support_exact(
        alternative = alt.u[i],
        probs = d,
        expectations = abs(seq_along(d) - 1 - lambda.u[i])
      )
      # pmin(1,
      #   switch(
      #     EXPR    = alternative,
      #     less    = c(cumsum(d[-length(d)]), 1),
      #     greater = c(1, rev(cumsum(rev(d[-1])))),
      #     minlike = ts.pv(statistics = d, probs = d),
      #     blaker  = ts.pv(
      #       statistics = pmin(
      #         c(cumsum(d[-length(d)]), 1),
      #         c(1, rev(cumsum(rev(d[-1]))))
      #       ),
      #       probs = d
      #     ),
      #     absdist = ts.pv(
      #       statistics = abs(seq_along(d) - 1 - lambda.u[i]),
      #       probs = d,
      #       decreasing = TRUE
      #     ),
      #     central = 2 * pmin(
      #       c(cumsum(d[-length(d)]), 1),
      #       c(1, rev(cumsum(rev(d[-1]))))
      #     )
      #   )
      # )
    } else {
      # find quantile with (approx.) smallest probability > 0
      q <- -floor(qnorm(2^-1074, -lambda.u[i], sqrt(lambda.u[i])))
      # greatest observation
      M <- max(q, N)
      # compute p-value support
      if(lambda.u[i] == 0){
        pv.supp <- switch(
          EXPR      = alt.u[i],
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
    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Exact Poisson test",
        paste0(
          "Normal-approximated Poisson test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(
        observations = data.frame(
          `number of events` = x,
          check.names = FALSE
        ),
        nullvalues = data.frame(
          `event rate` = lambda,
          check.names = FALSE
        ),
        parameters = data.frame(alternative = alternative)
      ),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = sapply(match.call(), deparse1)["x"]
    )
  } else res

  return(out)
}
