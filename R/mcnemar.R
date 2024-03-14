#' @title
#' McNemar's Test for Count Data
#'
#' @description
#' Performs McNemar's chi-square test or an exact variation of it to assess the
#' symmetry of rows and columns in a 2-by-2 contingency table. Only fourfold
#' tables can be processed, but multiple tables can be passed at once. The exact
#' version of the test is essentially an exact binomial test.
#'
#' @param x              an integer vector with four elements, a 2-by-2 matrix
#'                       or an integer matrix (or data frame) with four columns
#'                       where each line represents a 2-by-2 table to be tested,
#'                       i.e. a testing scenario.
#' @param alternative    indicates the alternative hypothesis and must be one of
#'                       \code{"two.sided"} (the default), \code{"less"} or
#'                       \code{"greater"}.
#' @param ts.method      indicate the two-sided p-value computation method (if
#'                       \code{alternative = two.sided}) and must be one of
#'                       \code{"minlike"} (the default), \code{"blaker"},
#'                       \code{"absdist"} or \code{"central"} (see details).
#'                       Ignored, if \code{exact = FALSE}.
#' @param exact          logical value that indicates whether p-values are to be
#'                       calculated by exact computation (\code{exact = TRUE};
#'                       the default) or by a continuous approximation.
#' @param correct        logical vaue that indicates if Edwards' continuity
#'                       correction is to be applied. Ignored, if
#'                       \code{exact = TRUE}.
#' @param simple.output  logical value that indicates whether the support sets,
#'                       i.e. all attainable p-values of each testing scenario
#'                       are to be returned, too (see below).
#'
#' @details
#' It can be shown that McNemar's test is a special case of the binomial test.
#' Therefore, its computations are handled by \code{\link{binom.test.pv}}.
#'
#' @return
#' If \code{supports = TRUE}, a list is returned:
#' \item{\code{$p.values}}{a vector of computed p-values as described above.}
#' \item{\code{$supports}}{a list of vectors, each containing all attainable
#'                         p-values of the respective scenario.}
#' Otherwise, only the vector of computed p-values is returned.
#'
#' @references
#' Agresti, A. (2002). \emph{Categorical data analysis}, 2nd ed. New York: John
#'   Wiley & Sons. pp. 350–354.
#'
#' @seealso
#' \code{\link[stats]{mcnemar.test}}, \code{\link{binom.test.pv}}
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
#' # Construction of exact p-values and their supports
#' results.ex  <- mcnemar.test.pv(df, simple.output = FALSE)
#' raw.pvalues <- results.ex$get_pvalues()
#' pCDFlist    <- results.ex$get_support_values()
#'
#' # Construction of chisquare p-values and their supports
#' #results.cs  <- mcnemar.test.pv(df, exact = FALSE, simple.output = FALSE)
#' #raw.pvalues <- results.cs$get_pvalues()
#' #pCDFlist    <- results.cs$get_support_values()
#'
#' @importFrom stats pchisq
#' @export
mcnemar.test.pv <- function(x, alternative = c("two.sided", "less", "greater"), ts.method = c("minlike", "blaker", "absdist", "central"), exact = TRUE, correct = TRUE, simple.output = TRUE){
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
  if(is.matrix(x)){
    # check if all values are finite, positive and integer
    if (any(!is.finite(x) | x < 0 | abs(x - (xr <- round(x))) > 1e-14))
      stop("All values of 'x' must be finite, non-negative and integer!")
    # round to integer
    x <- xr
    # stop immediately, if dimensions are violated
    if(all(dim(x) != c(2, 2)) && ncol(x) != 4 && nrow(x) != 4) stop(error.msg.x)
    # 2-by-2 matrices are transformed to single-row matrix
    if(all(dim(x) == c(2, 2))) x <- matrix(as.vector(x), 1, 4) else
      # transpose 4-row matrix (with more or less columns than 4) to 4-column matrix
      if((nrow(x) == 4 && ncol(x) != 4)) x <- t(x)
  }

  b <- x[, 2]
  c <- x[, 3]
#  n <- b + c

#  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

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
  return(binom.test.pv(b, b + c, 0.5, alternative, ts.method, exact, correct, simple.output))
}
