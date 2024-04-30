#' @title
#' Discrete Test Results Summary Class
#'
#' @description
#' This is the class used by `DiscreteTests` for summarising
#' [`DiscreteTestResults`] objects. It contains the summarised objects itself, as
#' well as a summary data frame as private members. Both can be read by public
#' methods.
#'
#' @examples
#' obj <- binom.test.pv(0:5, 5, 0.5)
#' res <- DiscreteTestResultsSummary$new(obj)
#' print(res)
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_class
#' @export
DiscreteTestResultsSummary <- R6Class(
  "summary.DiscreteTestResults",

  ## public ----

  public = list(
    #' @description
    #' Creates a new `summary.DiscreteTestResults` object.
    #'
    #' @param test_results   the [`DiscreteTestResults`] class object to be
    #'                       summarised.
    initialize = function(test_results) {
      # make sure that the results object is of class 'DiscreteTestResults'
      assert_class(test_results, c("DiscreteTestResults", "R6"))

      # create summary table
      inputs <- test_results$get_inputs(unique = FALSE)

      # compile data as list
      summary_table <- list(
        inputs$observations,
        inputs$nullvalues,
        inputs$parameters,
        test_results$get_pvalues()
      )

      # remove nulls and make data.frame
      summary_table <- as.data.frame(
        summary_table[!sapply(summary_table, is.null)]
      )

      # set column headers
      names(summary_table) <- c(
        names(inputs$observations),
        names(inputs$nullvalues),
        names(inputs$parameters),
        "p-value"
      )

      # assign inputs
      private$test_results      <- test_results
      private$summary_table     <- summary_table
    },

    #' @description
    #' Returns the underlying [DiscreteTestResults] object.
    #'
    #' @return
    #' A [DiscreteTestResults] R6 class object.
    get_test_results = function(){
      return(private$test_results)
    },

    #' @description
    #' Returns the summary table of the underlying [DiscreteTestResults] object.
    #' @return
    #' A data frame.
    get_summary_table = function(){
      return(private$summary_table)
    },

    #' @description
    #' Prints the summary.
    #'
    #' @param ...  further arguments passed to `print.data.frame`.
    #'
    #' @return
    #' Prints a summary table of the tested null hypotheses. The object itself
    #' is invisibly returned.
    print = function(...){
      print(private$test_results, FALSE, FALSE, FALSE)
      print(private$summary_table, ...)
      cat("\n")

      self
    }
  ),

  ## private ----

  private = list(
    # DiscreteTestResults R6 class object
    test_results = NULL,

    # summary table (a data frame)
    summary_table = data.frame()
  )
)

#' @title
#' Summarizing Discrete Test Results
#'
#' @description
#' `summary` method for class [`DiscreteTestResults`].
#'
#' @param object   object of class [`DiscreteTestResults`] to be summarised;
#'                 usually created by using one of the packages test functions,
#'                 e.g. [binom.test.pv()], with `simple.output = FALSE`.
#' @param ...      further arguments passed to or from other methods.
#'
#' @return
#' A [`summary.DiscreteTestResults`][DiscreteTestResultsSummary] R6 class
#' object.
#'
#' @examples
#' obj <- binom.test.pv(0:5, 5, 0.5)
#' summary(obj)
#'
#' @importFrom checkmate assert_class
#' @export
## S3 method for class 'DiscreteTestResults'
summary.DiscreteTestResults <- function(object, ...){
  assert_class(object, c("DiscreteTestResults", "R6"))

  DiscreteTestResultsSummary$new(object)
}
