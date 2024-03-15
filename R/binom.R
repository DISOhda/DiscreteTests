#' @title
#' Binomial Tests
#'
#' @description
#' `binom.test.pv` performs an exact binomial test or a normal
#' approximation about the probability of success in a Bernoulli experiment. It
#' is a vectorised version of `binom.test` that calculates only p-values.
#' Multiple testing scenarios can be passed at once. For two-sided ones, various
#' exact p-value methods are available.
#'
#' @param x              an integer vector giving the number of successes.
#' @param n              integer vector giving the number of trials.
#' @param p              hypothesized probabilities of success.
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x`, `n` and `p` are vectorised. They are
#' replicated automatically to have the same lengths. This allows multiple
#' hypotheses under the same conditions to be specified simultaneously.
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
#' # Construction of exact two-sided p-values ("blaker") and their supports
#' results.ex  <- binom.test.pv(k, n, p, ts.method = "blaker")
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_scenario_supports()
#'
#' # Construction of normal-approximated one-sided p-values (less) and their supports
#' results.ap  <- binom.test.pv(k, n, p, "less", exact = FALSE)
#' raw.pvalues <- results.ap$get_pvalues()
#' pCDFlist    <- results.ap$get_scenario_supports()
#'
#' @importFrom stats pnorm
#' @export
binom.test.pv <- function(x, n, p = 0.5, alternative = c("two.sided", "less", "greater"), ts.method = c("minlike", "blaker", "absdist", "central"), exact = TRUE, correct = TRUE, simple.output = FALSE){
  # plausability checks of input parameters
  len.x <- length(x)
  len.n <- length(n)
  len.p <- length(p)
  len.g <- max(len.x, len.n, len.p)

  if(is.list(x) || !is.vector(x) || !is.numeric(x))
    stop("'x' must be a numeric vector")

  if(min(len.x, len.n) < 1L)
    stop("Not enough data!")

  if(is.null(p)) p <- 0.5

  if(!is.null(p) && (is.list(p) || !is.vector(p) || !is.numeric(p) || any(is.na(p) | p < 0 | p > 1)))
    stop("'p' must be a vector of real numbers between 0 and 1")

  if(len.p < len.g) p <- rep_len(p, len.g)

  if(is.list(n) || !is.vector(n) || !is.numeric(n) || any(is.na(n) | !is.finite(n) | n < 0 | abs(n - (nr <- round(n))) > 1e-14))
    stop("'n' must be a vector of finite and non-negative integers")
  n <- nr

  if(len.n < len.g) n <- rep_len(n, len.g)

  if (any(is.na(x) | x < 0 | x > n | abs(x - (xr <- round(x))) > 1e-14))
    stop("All values of 'x' have to be non-negative, integer and must not exceed 'n'")
  x <- xr
  if(len.x < len.g) x <- rep_len(x, len.g)

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if(exact && alternative == "two.sided")
    alternative <- match.arg(ts.method, c("minlike", "blaker", "absdist", "central"))

  params <- unique(data.frame(n, p))
  len.u <- nrow(params)
  n.u <- params$n
  p.u <- params$p

  res <- numeric(len.g)
  if(!simple.output){
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  for(i in 1:len.u){
    idx <- which(n == n.u[i] & p == p.u[i])

    if(exact){
      d <- generate.binom.probs(n.u[i], p.u[i])
      pv.supp <- pmin(1, switch(alternative,
                                less    = cumsum(d),
                                greater = rev(cumsum(rev(d))),
                                minlike = ts.pv(d, d),
                                blaker  = ts.pv(pmin(cumsum(d), rev(cumsum(rev(d)))), d),
                                absdist = ts.pv(abs(0:n.u[i] - n.u[i] * p.u[i]), d, decreasing = TRUE),
                                central = 2 * pmin(cumsum(d), rev(cumsum(rev(d))))))
    }else{
      if(p.u[i] == 0)
        pv.supp <- switch(alternative,
                          less      = c(rep(1, n.u[i] + 1)),
                          greater   = c(1, rep(0, n.u[i])),
                          two.sided = c(1, rep(0, n.u[i])))
      else if(p.u[i] == 1)
        pv.supp <- switch(alternative,
                          less      = c(rep(0, n.u[i]), 1),
                          greater   = rep(1, n.u[i] + 1),
                          two.sided = c(rep(0, n.u[i]), 1))
      else{
        z <- 0:n.u[i] - n.u[i] * p.u[i]
        std <- sqrt(n.u[i] * p.u[i] * (1 - p.u[i]))
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
      ifelse(exact, "Exact binomial test", paste0("Normal-approximated binomial test", ifelse(correct, " with continuity correction", ""))),
      list(
        `number of successes` = x,
        Parameters = list(
          `number of trials` = n.u,
          `tested probability of success` = p.u
        )
      ),
      alternative,
      res,
      supports,
      indices,
      dname
    )
  }else res

  return(out)
}
