#' @title
#' McNemar's Test for Count Data
#'
#' @description
#' Performs McNemar's chi-square test or an exact variant to assess the symmetry
#' of rows and columns in a 2-by-2 contingency table. In contrast to
#' [stats::mcnemar.test()], it is vectorised, only calculates p-values and
#' offers their exact computation. Furthermore, it is capable of returning the
#' discrete p-value supports, i.e. all observable p-values under a null
#' hypothesis. Multiple tables can be analysed simultaneously. In two-sided
#' tests, several procedures of obtaining the respective p-values are
#' implemented. It is a special case of the [binomial test][binom.test.pv()].
#'
#' @param x   integer vector with four elements, a 2-by-2 matrix or an integer
#'            matrix (or data frame) with four columns where each line
#'            represents a 2-by-2 table to be tested.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts.method FALSE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' The parameters `x` and `alternative` are vectorised. They are replicated
#' automatically, such that the number of `x`'s rows is the same as the length
#' of `alternative`. This allows multiple hypotheses to be tested
#' simultaneously. Since `x is (if necessary) coerced to a matrix with four
#' columns, it is replicated row-wise.
#'
#' It can be shown that McNemar's test is a special case of the binomial test.
#' Therefore, its computations are handled by [binom.test.pv()]. In contrast to
#' that function, `mcnemar.test.pv` does not allow specifying exact two-sided
#' p-value calculation procedures. The reason is that McNemar's exact test
#' always tests for a probability of 0.5, in which case all these exact
#' two-sided p-value computation methods yield exactly the same results.
#'
#' @template return
#'
#' @references
#' Agresti, A. (2002). *Categorical data analysis* (2nd ed.). New York: John
#'   Wiley & Sons. pp. 411–413. \doi{10.1002/0471249688}
#'
#' @seealso
#' [stats::mcnemar.test()], [binom.test.pv()]
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
#' # Computation of exact p-values and their supports
#' results.ex  <- mcnemar.test.pv(df)
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_pvalue_supports()
#'
#' # Computation of chisquare p-values and their supports
#' results.cs  <- mcnemar.test.pv(df, exact = FALSE)
#' raw.pvalues <- results.cs$get_pvalues()
#' pCDFlist    <- results.cs$get_pvalue_supports()
#'
#' @importFrom stats pchisq
#' @importFrom checkmate assert_integerish
#' @export
mcnemar.test.pv <- function(
  x,
  alternative = "two.sided",
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
          make.names(paste(rep(colnames(x), rep(2,2)), rownames(x)))
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

  for(i in seq_len(len.a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
  }
  if(len.a < len.g) alternative <- rep_len(alternative, len.g)

  if(len.x < len.g)
    x <- matrix(rep_len(as.vector(t(x)), 4 * len.g), len.g, 4, TRUE)

  qassert(simple.output, "B1")

  b      <- x[, 2]
  c      <- x[, 3]
  n      <- b + c
  params <- unique(data.frame(n, alternative))
  len.u  <- nrow(params)
  alt.u  <- params$alternative

#  b <- x[, 2]
##  c <- x[, 3]
#  n <- b + x[, 3]
##  n <- b + c

#  if(exact || (!exact && alternative != "two.sided")){
#    res <- binom.test.pv(a, n, 0.5, alternative, "central", exact, correct, supports)
#  }else{
#    D <- b - c
#    #chi <- (abs(D) - (D & correct))^2 / n
#    pv <- switch(alternative,
#                 less = pnorm(D, -correct, sqrt(n)), #pnorm(D, sign(D) * correct, sqrt(n)),
#                 greater = pnorm(D, correct, sqrt(n), lower.tail = FALSE), #pnorm(D, sign(D) * correct, sqrt(n), lower.tail = FALSE),
#                 two.sided = 2 * pnorm(-abs(D), -correct, sqrt(n)))
#                 #two.sided = pchisq((abs(D) - (D & correct))^2 / n, 1, lower.tail = FALSE))
#    zeros <- which(n == 0)
#    if(length(zeros)) pv[zeros] <- 1
#
#    if(simple.output){
#      res <- pv
#    }else{
#      num <- length(n)
#      res <- list(p.values = NULL, supports = NULL)
#      res$p.values <- numeric(num)
#      res$supports <- vector("list", num)
#
#      n.u <- unique(n)
#      len <- length(n.u)
#      for(i in 1:len){
#        if(n.u[i]){
#          D <- seq(-n.u[i], n.u[i], 2)
#          #chi <- (abs(D) - (D & correct))^2 / n.u[i]
#          supp <- switch(alternative,
#                         less = pnorm(D, -correct, sqrt(n.u[i])), #pnorm(D, sign(D) * correct, sqrt(n.u[i])),
#                         greater = pnorm(D, correct, sqrt(n.u[i]), lower.tail = FALSE), #pnorm(D, sign(D) * correct, sqrt(n.u[i]), lower.tail = FALSE),
#                         two.sided = 2 * pnorm(-abs(D), -correct, 0, sqrt(n.u[i])))
#                         #two.sided = pchisq((abs(D) - (D & correct))^2 / n.u[i], 1, lower.tail = FALSE))
#        }else supp <- 1#
#
#        idx <- which(n == n.u[i])
#        for(j in idx){
#          res$supports[[j]] <- sort(unique(supp))
#        }
#      }
#      res$p.values <- pv
#    }
#  }
#
#  return(res)

  res <- binom.test.pv(b, n, 0.5, alternative, "central", exact, correct, simple.output)

  out <- if(!simple.output) {
    if(is.null(colnames(x)))
      colnames(x) <- paste0("x", c("[1, 1]", "[2, 1]", "[1, 2]", "[2, 2]"))

    pars <- res$get_inputs(unique = TRUE)

    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "McNemar's exact test",
        paste0(
          "McNemar's Chi-squared Test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = c(
        list(
          observations = as.data.frame(x),
          nullvalues = data.frame(
            `counter-diagonal values proportions` = rep(0.5, len.g),
            check.names = FALSE
          ),
          parameters = data.frame(
            `counter-diagonal sum` = n,
            alternative = alternative,
            check.names = FALSE
          )
        )
      ),
      p_values = res$get_pvalues(),
      pvalue_supports = res$get_pvalue_supports(unique = TRUE),
      support_indices = res$get_support_indices(),
      data_name = sapply(match.call(), deparse1)["x"]
    )
  } else res

  return(out)
}
