#' @name fisher_test_pv
#'
#' @title
#' Fisher's Exact Test for Count Data
#'
#' @description
#' `fisher_test_pv()` performs Fisher's exact test or a chi-square approximation
#' to assess if rows and columns of a 2-by-2 contingency table with fixed
#' marginals are independent. In contrast to [`stats::fisher.test()`], it is
#' vectorised, only calculates *p*-values and offers a normal approximation of
#' their computation. Furthermore, it is capable of returning the discrete
#' *p*-value supports, i.e. all observable *p*-values under a null hypothesis.
#' Multiple tables can be analysed simultaneously. In two-sided tests, several
#' procedures of obtaining the respective *p*-values are implemented.
#'
#' `r lifecycle::badge('deprecated')`\cr
#' **Note**: Please use `fisher_test_pv()`! The older `fisher.test.pv()` is
#' deprecated in order to migrate to snake case. It will be removed in a future
#' version.
#'
#' @param x   integer vector with four elements, a 2-by-2 matrix or an integer
#'            matrix (or data frame) with four columns, where each line
#'            represents a 2-by-2 table to be tested.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts_method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x` and `alternative` are vectorised. They are replicated
#' automatically, such that the number of `x`'s rows is the same as the length
#' of `alternative`. This allows multiple null hypotheses to be tested
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
#' [`stats::fisher.test()`]
#'
#' @references
#' Fisher, R. A. (1935). The logic of inductive inference.
#'   *Journal of the Royal Statistical Society Series A*, **98**, pp.
#'   39–54. \doi{10.2307/2342435}
#'
#' Agresti, A. (2002). *Categorical data analysis*. Second Edition. New York:
#'   John Wiley & Sons. pp. 91–97. \doi{10.1002/0471249688}
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
#' S1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' S2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#' N1 <- rep(148, 9)
#' N2 <- rep(132, 9)
#' F1 <- N1 - S1
#' F2 <- N2 - S2
#' df <- data.frame(S1, F1, S2, F2)
#'
#' # Fisher's exact p-values and their supports
#' results_f   <- fisher_test_pv(df)
#' raw_pvalues <- results_f$get_pvalues()
#' pCDFlist    <- results_f$get_pvalue_supports()
#'
#' # Chi-square-approximated p-values and their supports
#' results_c   <- fisher_test_pv(df, exact = FALSE)
#' raw_pvalues <- results_c$get_pvalues()
#' pCDFlist    <- results_c$get_pvalue_supports()
#'
#' @importFrom stats dhyper pchisq
#' @importFrom checkmate assert_integerish qassert
#' @export
fisher_test_pv <- function(
  x,
  alternative = "two.sided",
  ts_method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # plausibility checks of input parameters

  # define error message for malformed x
  error_msg_x <- paste("'x' must either be a 2-by-2 matrix,",
                       "a four-element vector or a four-column matrix")

  #  if x is a vector, make it a matrix with one row
  if(is.vector(x) && !is.list(x))
    x <- t(x)
  # if x is a data frame, make it a matrix
  if(is.data.frame(x))
    x <- as.matrix(x)
  # if x is a list, then abort
  if(is.list(x)) stop(error_msg_x)
  # when x is a matrix, it must satisfy some conditions
  if(is.matrix(x)) {
    # check if all values are non-negative and close to integer
    assert_integerish(x, lower = 0)
    # round to integer
    x <- round(x)
    # stop immediately, if dimensions are violated
    if(any(dim(x) != c(2, 2)) && ncol(x) != 4 && nrow(x) != 4)
      stop(error_msg_x)
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
  } else stop(error_msg_x)
  len_x <- nrow(x)

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

  # determine parameters for R's hypergeometric distribution implementation
  m <- x[, 1] + x[, 2] # sums of 1st columns
  n <- x[, 3] + x[, 4] # sums of 2nd columns
  k <- x[, 1] + x[, 3] # sums of 1st rows
  q <- x[, 1]          # upper left elements

  # replicate inputs to same length
  len_g <- max(len_x, len_a)
  if(len_x < len_g){
    m <- rep_len(m, len_g)
    n <- rep_len(n, len_g)
    k <- rep_len(k, len_g)
    q <- rep_len(q, len_g)
  }
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # determine unique parameter sets and possible "q" value boundaries (support)
  params <- unique(data.frame(m, n, k, alternative))
  len_u <- nrow(params)
  m_u <- as.integer(params$m)
  n_u <- as.integer(params$n)
  k_u <- as.integer(params$k)
  alt_u <- params$alternative
  hi <- pmin(k_u, m_u)
  lo <- pmax(0, k_u - n_u)

  # prepare output
  res <- numeric(nrow(x))
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # begin computations
  for(i in seq_len(len_u)) {
    # which hypotheses belong to the current unique parameter set
    idx_supp <- which(m == m_u[i] & n == n_u[i] & k == k_u[i] & alternative == alt_u[i])
    # possible "q" values
    support <- lo[i]:hi[i]

    if(exact){
      # compute all possible probability masses under fixed marginals
      d <- numerical_adjust(dhyper(support, m_u[i], n_u[i], k_u[i]))
      # p-value supports according to alternative and (maybe) two-sided method
      pv_supp <- support_exact(
        alternative = alt_u[i],
        probs = d,
        expectations = abs(support - m_u[i] * k_u[i] / (n_u[i] + m_u[i]))
      )
    } else {
      # observable tables under fixed marginals
      y <- rbind(
        support,
        m_u[i] - support,
        k_u[i] - support,
        n_u[i] + support - k_u[i]
      )
      # observable odds ratios
      obs_or <- y[1, ] * y[4, ] / (y[2, ] * y[3, ])
      # direction of deviance from homogeneity (odds ratio = 1)
      delta <- sign(obs_or - 1)
      # correct for NaN deltas (if observable odds ratio = NaN)
      delta[is.nan(delta)] <- 0

      # expected table under homogeneity
      f11 <- m_u[i] * k_u[i] / (n_u[i] + m_u[i])
      f10 <- k_u[i] - f11
      f01 <- m_u[i] - f11
      f00 <- n_u[i] - f10
      expected <- pmax(0, c(f11, f01, f10, f00))
      if(any(expected < 5))
        warning("One or more Chi-squared approximations may be incorrect!\n")

      # chi-square values
      absdiff <- if(correct) {
        matrix(pmax(0, abs(y - expected) - 0.5), 4)^2
      } else (y - expected)^2
      chi <- absdiff / expected
      chi[is.nan(chi)] <- 0
      chi <- numerical_adjust(colSums(chi), FALSE)
      # degrees of freedom
      df <- 1 - any(expected == 0)
      # p-value supports according to alternative
      pv_supp <- switch(alt_u[i],
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

    # store results and support
    idx_obs <- sapply(seq_along(idx_supp), function(j) which(support == q[idx_supp[j]]))
    res[idx_supp] <- pv_supp[idx_obs]
    if(!simple_output) {
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)
    if(is.null(colnames(x)))
      colnames(x) <- paste0("x", c("[1, 1]", "[2, 1]", "[1, 2]", "[2, 2]"))
    if(len_x < len_g)
      x <- x[rep_len(seq_len(len_x), len_g), ]

    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Fisher's Exact test",
        "Chi-squared test for homogenity"
      ),
      inputs = list(
        observations = as.data.frame(x),
        nullvalues = data.frame(
          `odds ratio` = rep(1, len_g),
          check.names = FALSE
        ),
        parameters = data.frame(
          `first column sum` = m,
          `second column sum` = n,
          `first row sum` = k,
          alternative = alternative,
          exact = exact,
          distribution = ifelse(exact, "hypergeometric",
            ifelse(alternative == "two.sided", "chi-squared", "normal")
          ),
          check.names = FALSE
        )
      ),
      statistics = NULL,
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = dnames["x"]
    )
  } else res

  return(out)
}

#' @rdname fisher_test_pv
#' @export
#' @importFrom lifecycle deprecate_soft
fisher.test.pv <- function(
    x,
    alternative = "two.sided",
    ts.method = "minlike",
    exact = TRUE,
    correct = TRUE,
    simple.output = FALSE
) {
  deprecate_soft("0.2", "fisher.test.pv()", "fisher_test_pv()")

  fisher_test_pv(x, alternative, ts.method, exact, correct, simple.output)
}
