#' @title
#' Fisher's Exact Test for Count Data
#'
#' @description
#' `fisher.test.pv()` performs Fisher's exact test or a chi-square approximation
#' to assess if rows and columns of a 2-by-2 contingency table with fixed
#' marginals are independent. In contrast to [stats::fisher.test()], it is
#' vectorised, only calculates p-values and offers a normal approximation of
#' their computation. Furthermore, it is capable of returning the discrete
#' p-value supports, i.e. all observable p-values under a null hypothesis.
#' Multiple tables can be analysed simultaneously. In two-sided tests, several
#' procedures of obtaining the respective p-values are implemented.
#'
#' @param x   integer vector with four elements, a 2-by-2 matrix or an integer
#'            matrix (or data frame) with four columns, where each line
#'            represents a 2-by-2 table to be tested.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x` and `alternative` are vectorised. They are replicated
#' automatically, such that the number of `x`'s rows is the same as the length
#' of `alternative`. This allows multiple hypotheses to be tested
#' simultaneously. Since `x` is (if necessary) coerced to a matrix with four
#' columns, it is replicated row-wise.
#'
#' If `exact = TRUE`, Fisher's exact test is performed (the specific hypothesis
#' depends on the value of `alternative`). Otherwise, if `exact = FALSE`, a
#' chi-square approximation is used for two-sided hypotheses or a normal
#' approximation for one-sided tests, based on the square root of the
#' chi-squared statistic. This is possible because the degrees of freedom of
#' chi-squared tests on 2-by-2 tables are limited to 1.
#'
#' @template details_two_sided
#'
#' @template return
#'
#' @seealso
#' [stats::fisher.test()]
#'
#' @references
#' Fisher, R. A. (1935). The logic of inductive inference.
#'   *Journal of the Royal Statistical Society Series A*, **98**, pp.
#'   39–54. \doi{10.2307/2342435}
#'
#' Agresti, A. (2002). *Categorical data analysis* (2nd ed.). New York: John
#'   Wiley & Sons. pp. 91–101.
#'
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. *Canadian Journal of Statistics*,
#'   **28**(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Hirji, K. F. (2006). *Exact analysis of discrete data*. New York: Chapman
#'   and Hall/CRC. pp. 55-83
#'
#' @examples
#' # Constructing
#' S1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' S2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#' N1 <- rep(148, 9)
#' N2 <- rep(132, 9)
#' F1 <- N1 - S1
#' F2 <- N2 - S2
#' df <- data.frame(S1, F1, S2, F2)
#'
#' # Computation of Fisher's exact p-values (default: "minlike") and their supports
#' results.f   <- fisher.test.pv(df)
#' raw.pvalues <- results.f$get_pvalues()
#' pCDFlist    <- results.f$get_pvalue_supports()
#'
#' # Computation of p-values of chi-square tests and their supports
#' results.c   <- fisher.test.pv(df, exact = FALSE)
#' raw.pvalues <- results.c$get_pvalues()
#' pCDFlist    <- results.c$get_pvalue_supports()
#'
#' @importFrom stats dhyper pnorm pchisq
#' @importFrom checkmate assert_integerish
#' @export
fisher.test.pv <- function(
  x,
  alternative = "two.sided",
  ts.method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple.output = FALSE
) {
  # plausibility checks of input parameters

  # define error message for malformed x
  error.msg.x <- paste("'x' must either be a 2-by-2 matrix,",
                       "a four-element vector or a four-column matrix")

  #  if x is a vector, make it a matrix with one row
  if(is.vector(x) && !is.list(x))
    x <- t(x)
  # if x is a data frame, make it a matrix
  if(is.data.frame(x))
    x <- as.matrix(x)
  # if x is a list, then abort
  if(is.list(x)) stop(error.msg.x)
  # when x is a matrix, it must satisfy some conditions
  if(is.matrix(x)) {
    # check if all values are non-negative and close to integer
    assert_integerish(x, lower = 0)
    # round to integer
    x <- round(x)
    # stop immediately, if dimensions are violated
    if(all(dim(x) != c(2, 2)) && ncol(x) != 4 && nrow(x) != 4)
      stop(error.msg.x)
    # 2-by-2 matrices are transformed to single-row matrix
    if(all(dim(x) == c(2, 2))) {
      x <- matrix(as.vector(x), 1, 4,
        dimnames = list(NULL,
          make.names(paste(rep(colnames(x), rep(2, 2)), rownames(x)))
        )
      )
    } else
      # transpose 4-row matrix (with more or less columns than 4) to 4-column matrix
      if((nrow(x) == 4 && ncol(x) != 4))
        x <- t(x)
  } else stop(error.msg.x)

  # lengths
  len.x <- nrow(x)
  len.a <- length(alternative)
  len.g <- max(len.x, len.a)

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

  ## computations
  #  parameters for R's hypergeometric distribution implementation
  m <- x[, 1] + x[, 2] # sums of 1st columns
  n <- x[, 3] + x[, 4] # sums of 2nd columns
  k <- x[, 1] + x[, 3] # sums of 1st rows
  q <- x[, 1]          # upper left elements

  # recycle, if necessary
  if(len.x < len.g){
    m <- rep_len(m, len.g)
    n <- rep_len(n, len.g)
    k <- rep_len(k, len.g)
    q <- rep_len(q, len.g)
  }

  # determine unique parameter sets and possible "q" value boundaries (support)
  params <- unique(data.frame(m, n, k, alternative))
  len.u <- nrow(params)
  m.u <- as.integer(params$m)
  n.u <- as.integer(params$n)
  k.u <- as.integer(params$k)
  alt.u <- params$alternative
  hi <- pmin(k.u, m.u)
  lo <- pmax(0, k.u - n.u)

  # prepare output
  res <- numeric(nrow(x))
  if(!simple.output) {
    supports <- vector("list", len.u)
    indices  <- vector("list", len.u)
  }

  # loop through unique parameter sets
  for(i in 1:len.u) {
    # which hypotheses belong to the current unique parameter set
    idx <- which(m == m.u[i] & n == n.u[i] & k == k.u[i] & alternative == alt.u[i])
    # possible "q" values
    support <- lo[i]:hi[i]

    if(exact){
      # compute all possible probability masses under fixed marginals
      d <- numerical.adjust(dhyper(support, m.u[i], n.u[i], k.u[i]))
      # p-value supports according to alternative and (maybe) two-sided method
      pv.supp <- support_exact(
        alternative = alt.u[i],
        probs = d,
        expectations = abs(support - m.u[i] * k.u[i] / (n.u[i] + m.u[i]))
      )
      # pv.supp <- pmin(1,
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
      #       statistics = abs(support - m.u[i] * k.u[i] / (n.u[i] + m.u[i])),
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
      # observable tables under fixed marginals
      y <- rbind(
        support,
        m.u[i] - support,
        k.u[i] - support,
        n.u[i] + support - k.u[i]
      )
      # observable odds ratios
      obs.or <- y[1, ] * y[4, ] / (y[2, ] * y[3, ])
      # direction of deviance from homogeneity (odds ratio = 1)
      delta <- sign(obs.or - 1)
      # correct for NaN deltas (if observable odds ratio = NaN)
      delta[is.nan(delta)] <- 0

      # expected table under homogeneity
      f11 <- m.u[i] * k.u[i] / (n.u[i] + m.u[i])
      f10 <- k.u[i] - f11
      f01 <- m.u[i] - f11
      f00 <- n.u[i] - f10
      expected <- pmax(0, c(f11, f01, f10, f00))
      if(any(expected < 5))
        warning("One or more Chi-squared approximations may be incorrect!\n")

      # chi-square values
      absdiff <- if(correct) {
        matrix(pmax(0, abs(y - expected) - 0.5), 4)^2
      } else (y - expected)^2
      chi <- absdiff / expected
      chi[is.nan(chi)] <- 0
      chi <- numerical.adjust(colSums(chi), FALSE)
      # degrees of freedom
      df <- 1 - any(expected == 0)
      # p-value supports according to alternative
      pv.supp <- switch(alt.u[i],
        less = {
          p <- pnorm_zero(delta * sqrt(chi), df)
          p[p == max(p)] <- 1
          p
        },
        greater = {
          p <- pnorm_zero(delta * sqrt(chi), df, FALSE)
          p[p == max(p)] <- 1
          p
        },
        two.sided = pchisq(chi, df, lower.tail = FALSE)
      )
    }

    idx.obs <- sapply(seq_along(idx), function(j) which(support == q[idx[j]]))
    res[idx] <- pv.supp[idx.obs]
    if(!simple.output) {
      supports[[i]] <- sort(unique(pv.supp))
      indices[[i]]  <- idx
    }
  }

  out <- if(!simple.output) {
    if(is.null(colnames(x)))
      colnames(x) <- paste0("x", c("[1, 1]", "[2, 1]", "[1, 2]", "[2, 2]"))
    if(len.x < len.g)
      x <- x[rep_len(seq_len(len.x), len.g), ]

    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Fisher's Exact Test",
        paste0(
          "Chi-squared test for homogenity",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(
        observations = as.data.frame(x),
        nullvalues = data.frame(
          `odds ratio` = rep(1, len.g),
          check.names = FALSE
        ),
        parameters = data.frame(
          `first column sum` = m,
          `second column sum` = n,
          `first row sum` = k,
          alternative = alternative,
          check.names = FALSE
        )
      ),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = sapply(match.call(), deparse1)["x"]
    )
  } else res

  return(out)
}
