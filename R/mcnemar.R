#' @name mcnemar_test_pv
#'
#' @title
#' McNemar's Test for Count Data
#'
#' @description
#' Performs McNemar's chi-square test or an exact variant to assess the symmetry
#' of rows and columns in a 2-by-2 contingency table. In contrast to
#' [`stats::mcnemar.test()`], it is vectorised, only calculates *p*-values and
#' offers their exact computation. Furthermore, it is capable of returning the
#' discrete *p*-value supports, i.e. all observable *p*-values under a null
#' hypothesis. Multiple tables can be analysed simultaneously. In two-sided
#' tests, several procedures of obtaining the respective *p*-values are
#' implemented. It is a special case of the [binomial test][binom_test_pv()].
#'
#' `r lifecycle::badge('deprecated')`\cr
#' **Note**: Please use `mcnemar_test_pv()`! The older `mcnemar.test.pv()` is
#' deprecated in order to migrate to snake case. It will be removed in a future
#' version.
#'
#' @param x   integer vector with four elements, a 2-by-2 matrix or an integer
#'            matrix (or data frame) with four columns where each line
#'            represents a 2-by-2 table to be tested.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x` and `alternative` are vectorised. They are replicated
#' automatically, such that the number of `x`'s rows is the same as the length
#' of `alternative`. This allows multiple null hypotheses to be tested
#' simultaneously. Since `x is (if necessary) coerced to a matrix with four
#' columns, it is replicated row-wise.
#'
#' It can be shown that McNemar's test is a special case of the binomial test.
#' Therefore, its computations are handled by [`binom_test_pv()`]. In
#' contrast to that function, `mcnemar_test_pv()` does not allow specifying
#' exact two-sided p-value calculation procedures. The reason is that McNemar's
#' exact test always tests for a probability of 0.5, in which case all these
#' exact two-sided p-value computation methods yield exactly the same results.
#'
#' @template return
#'
#' @references
#' Agresti, A. (2002). *Categorical data analysis*. Second Edition. New York:
#'   John Wiley & Sons. pp. 411–413. \doi{10.1002/0471249688}
#'
#' @seealso
#' [`stats::mcnemar.test()`], [`binom_test_pv()`]
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
#' # Exact p-values and their supports
#' results_ex  <- mcnemar_test_pv(df)
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Chi-square-approximated p-values and their supports
#' results_cs  <- mcnemar_test_pv(df, exact = FALSE)
#' raw_pvalues <- results_cs$get_pvalues()
#' pCDFlist    <- results_cs$get_pvalue_supports()
#'
#' @importFrom checkmate assert_integerish qassert
#' @importFrom cli cli_abort
#' @export
mcnemar_test_pv <- function(
  x,
  alternative = "two.sided",
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
  if(is.list(x)) cli_abort(error_msg_x)
  # when x is a matrix, it must satisfy some conditions
  if(is.matrix(x)) {
    # check if all values are non-negative and close to integer
    assert_integerish(x, lower = 0)
    # round to integer
    x <- round(x)
    mode(x) <- "integer"
    # stop immediately, if dimensions are violated
    if(any(dim(x) != c(2, 2)) && ncol(x) != 4 && nrow(x) != 4)
      cli_abort(error_msg_x)
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
  } else cli_abort(error_msg_x)

  # lengths
  len_x <- nrow(x)
  len_a <- length(alternative)
  len_g <- max(len_x, len_a)

  # check exactness and continuity correction arguments
  qassert(exact, "B1")
  if(!exact) qassert(correct, "B1")

  # check alternatives
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
  }
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # enlarge input matrix if necessary
  if(len_x < len_g)
    x <- matrix(rep_len(as.vector(t(x)), 4 * len_g), len_g, 4, TRUE)

  # check output type argument
  qassert(simple_output, "B1")

  # test parameters
  b <- x[, 2]
  c <- x[, 3]
  n <- b + c

  # compute test results
  res <- binom_test_pv(b, n, 0.5, alternative, "central", exact, correct, simple_output)

  # create output object
  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)
    if(is.null(colnames(x)))
      colnames(x) <- paste0("x", c("[1, 1]", "[2, 1]", "[1, 2]", "[2, 2]"))

    DiscreteTestResults$new(
      test_name = "McNemar's test",
      inputs = c(
        list(
          observations = as.data.frame(x),
          nullvalues = data.frame(
            `counter-diagonal values proportion` = rep(0.5, len_g),
            check.names = FALSE
          ),
          parameters = NULL,
          computation = data.frame(
            alternative = alternative,
            exact = exact,
            distribution = ifelse(exact, "binomial", "normal"),
            distribution.size = if(exact) n else NA_integer_,
            distribution.mean = if(exact) NA_real_ else n * 0.5,
            distribution.sd = if(exact) NA_real_ else sqrt(n * 0.25),
            `continuity correction` = if(exact) NA else correct,
            check.names = FALSE
          )
        )
      ),
      statistics = NULL,
      p_values = res$get_pvalues(),
      pvalue_supports = res$get_pvalue_supports(unique = TRUE),
      support_indices = res$get_support_indices(),
      data_name = dnames["x"]
    )
  } else res

  # return results
  return(out)
}

#' @rdname mcnemar_test_pv
#' @export
#' @importFrom lifecycle deprecate_soft
mcnemar.test.pv <- function(
    x,
    alternative = "two.sided",
    exact = TRUE,
    correct = TRUE,
    simple.output = FALSE
) {
  deprecate_soft("0.2", "mcnemar.test.pv()", "mcnemar_test_pv()")

  mcnemar_test_pv(x, alternative, exact, correct, simple.output)
}
