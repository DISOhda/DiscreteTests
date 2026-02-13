#' @title Conditional Two-Sample Homogeneity Test for Binomial Experiments
#'
#' @description
#' Performs an exact or approximate conditional test about the homogeneity of
#' two binomial samples, i.e. regarding the respective probabilities of success.
#' It is vectorised, only calculates *p*-values and offers a normal
#' approximation of their computation. Furthermore, it is capable of returning
#' the discrete *p*-value supports, i.e. all observable *p*-values under a null
#' hypothesis. Multiple tests can be evaluated simultaneously. In two-sided
#' tests, several procedures of obtaining the respective *p*-values are
#' implemented.
#'
#' @param x   integer vector with two elements or a matrix with two columns or a
#'            data frame with two columns giving the number of successes for the
#'            two experiments.
#' @param n   integer vector with two elements or a matrix with two columns or a
#'            data frame with two columns giving the number of trials for the
#'            two experiments.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts_method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#'
#' @details
#' The parameters `x`, `n` and `alternative` are vectorised. They are replicated
#' automatically, such that the number of `x`'s rows is the same as the length
#' of `alternative`. This allows multiple null hypotheses to be tested
#' simultaneously. Since `x` and `n` are coerced to matrices (if necessary) with
#' two columns, they are replicated row-wise.
#'
#' It can be shown that this test is a special case of Fisher's exact test,
#' because it is conditional on the numbers of trials **and** the sums of
#' successes and failures. Therefore, its computations are handled by
#' [`fisher_test_pv()`].
#'
#' @template details_two_sided
#'
#' @template return
#'
#' @seealso
#' [`fisher_test_pv()`]
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
#' set.seed(3)
#' p1 <- c(0.25, 0.5, 0.75)
#' p2 <- c(0.15, 0.5, 0.60)
#' n1 <- c(10, 20, 50)
#' n2 <- c(25, 75, 200)
#' x1 <- rbinom(3, n1, p1)
#' x2 <- rbinom(3, n2, p2)
#' x  <- cbind(x1 = x1, x2 = x2)
#' n  <- cbind(n1 = n1, n2 = n2)
#'
#' # Exact two-sided p-values ("blaker") and their supports
#' results_ex <- homogeneity_test_pv(x, n, ts_method = "blaker")
#' print(results_ex)
#' results_ex$get_pvalues()
#' results_ex$get_pvalue_supports()
#'
#' # Normal-approximated one-sided p-values ("less") and their supports
#' results_ap <- homogeneity_test_pv(x, n, "less", exact = FALSE)
#' print(results_ap)
#' results_ap$get_pvalues()
#' results_ap$get_pvalue_supports()
#'
#' @importFrom checkmate assert_integerish
#' @importFrom cli cli_abort
#' @export
homogeneity_test_pv <- function(
  x,
  n,
  alternative = "two.sided",
  ts_method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # plausibility checks of input parameters

  # define error message for malformed x or n
  error_msg <- paste("must either be a two-element vector, a two-column matrix",
                     "or a two-column data frame")

  # make sure that x is a data frame, matrix or vector
  qassert(x, c("D+", "M+", "V2"))
  # make sure that n is a data frame, matrix or vector
  qassert(n, c("D+", "M+", "V2"))

  # if x is an atomic vector, make it a matrix with one row
  if(is.vector(x) && is.atomic(x))
    x <- t(x)
  # if n is an atomic vector, make it a matrix with one row
  if(is.vector(n) && is.atomic(n))
    n <- t(n)

  # if x is a data frame, make it a matrix
  if(is.data.frame(x))
    x <- as.matrix(x)
  # if n is a data frame, make it a matrix
  if(is.data.frame(n))
    n <- as.matrix(n)

  # when x is a matrix, it must satisfy some conditions
  if(is.matrix(x)) {
    # stop immediately, if dimensions are violated
    if(ncol(x) != 2 && nrow(x) != 2)
      cli_abort(paste("'x'", error_msg))
    # check if all values are non-negative and close to integer
    assert_integerish(x, lower = 0, min.len = 2)
    # coerce to integer
    x <- round(x)
    mode(x) <- "integer"
  } else cli_abort(paste("'x'", error_msg))
  # when n is a matrix, it must satisfy some conditions
  if(is.matrix(n)) {
    # stop immediately, if dimensions are violated
    if(ncol(n) != 2 && nrow(n) != 2)
      cli_abort(paste("'n'", error_msg))
    # check if all values are non-negative and close to integer
    assert_integerish(n, lower = 0, min.len = 2)
    # coerce to integer
    n <- round(n)
    mode(n) <- "integer"
  } else cli_abort(paste("'n'", error_msg))

  # determine largest number of rows
  len_x <- nrow(x)
  len_n <- nrow(n)
  len_g <- max(len_x, len_n)
  # enlarge x by row-replication, if necessary
  if(len_x < len_g)
    x <- x[rep_len(seq_len(len_x), len_g), ]
  # enlarge n by row-replication, if necessary
  if(len_n < len_g)
    n <- n[rep_len(seq_len(len_n), len_g), ]

  # create matrix of vectorised 2-by-2 tables
  y <- cbind(x, n - x);

  # compute test results
  res <- fisher_test_pv(y, alternative, ts_method, exact, correct, simple_output)

  # create output object
  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)
    #if(is.null(colnames(x)))
      colnames(x) <- paste(c("first", "second"), "number of successes")
    #if(is.null(colnames(n)))
      colnames(n) <- paste(c("first", "second"), "number of trials")

    DiscreteTestResults$new(
      test_name = "Conditional two-sample binomial homogeneity test",
      inputs = c(
        list(
          observations = data.frame(x, check.names = FALSE),
          parameters = data.frame(n, check.names = FALSE),
          nullvalues = data.frame(
            `odds ratio` = rep(1, len_g),
            check.names = FALSE
          ),
          computation = Filter(
            function(df) !all(is.na(df)),
            data.frame(
              alternative = res$get_inputs()$computation$alternative,
              exact = res$get_inputs()$computation$exact,
              distribution = res$get_inputs()$computation$distribution,
              #distribution.df = ifelse(
              #  !exact & alternative == "two.sided", df, NA_integer_
              #),
              #distribution.mean = ifelse(
              #  !exact & alternative != "two.sided", 0, NA_real_
              #),
              #distribution.sd = ifelse(
              #  !exact & alternative != "two.sided", 1, NA_real_
              #),
              `continuity correction` = if(exact) NA else correct,
              `total number of successes` = rowSums(x[, 1:2, drop = FALSE]),
              `total number of trials` = rowSums(n),
              check.names = FALSE
            )
          )
        )
      ),
      statistics = res$get_statistics(),
      p_values = res$get_pvalues(),
      pvalue_supports = res$get_pvalue_supports(unique = TRUE),
      support_indices = res$get_support_indices(),
      data_name = dnames[c("x", "n")]
    )
  } else res

  # return results
  return(out)
}
