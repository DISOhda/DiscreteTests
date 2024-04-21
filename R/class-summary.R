#' @title
#' Discrete Test Results Summary
#'
#' @description
#' This is the class used by `DiscreteTests` for summarising
#' `DiscreteTestResults` objects. It contains the summarised objects itself, as
#' well as a summary data frame as private members. Both can be read by public
#' methods. **Note**: The class generator is not exported to the namespace and
#' must be accessed via `DiscreteTests:::DiscreteTestResultsSummary`.
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_class assert_data_frame
DiscreteTestResultsSummary <- R6Class(
  "summary.DiscreteTestResults",

  ## public ----

  public = list(
    #' @description
    #' Creates a new `summary.DiscreteTestResults` object.
    #'
    #' @param test_results   a [DiscreteTestResults] class object.
    initialize = function(test_results) {
      # make sure that the results object is of class 'DiscreteTestResults'
      assert_class(test_results, c("DiscreteTestResults", "R6"))

      # create summary table
      inputs <- test_results$get_inputs(unique = FALSE)

      summary_table <- list(
        inputs[[1]],
        inputs[[2]],
        test_results$get_pvalues()
      )
      summary_table <- as.data.frame(
        summary_table[!sapply(summary_table, is.null)]
      )
      if(is.matrix(inputs[[1]]) || is.data.frame(inputs[[1]]))
        name_obs <- colnames(inputs[[1]]) else
          name_obs <- names(inputs[1])
      names(summary_table) <- c(name_obs, names(inputs[[2]]), "p-value")

      # assign inputs
      private$test_results      <- test_results
      private$summary_table     <- summary_table
    },

    #' @description
    #' Returns the underlying [DiscreteTestResults] object.
    #' @returns
    #' A [DiscreteTestResults] R6 class object.
    get_test_results = function(){
      return(private$test_results)
    },

    #' @description
    #' Returns the summary table of the underlying [DiscreteTestResults] object.
    #' @returns
    #' A data frame.
    get_summary_table = function(){
      return(private$summary_table)
    },

    #' @description
    #' Prints the summary.
    #' @param ...  further arguments passed to `print.data.frame`.
    #' @returns
    #' Prints a summary of the tested null hypotheses. The object itself is
    #' invisibly returned.
    print = function(...){
      print(private$test_results, FALSE, FALSE, FALSE)
      cat("\n")
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
    summary_table = data.frame(),

    # version of class definition
    class_version = "0.1.0"
  )
)

#' @title
#' Summarizing Discrete Test Results
#'
#' @description
#' `summary` method for class `"DiscreteTestResults"`.
#'
#' @param object  an object of class `"DiscreteTestResults"`, usually
#'                produced by a call to one of the packages test functions, e.g.
#'                [binom.test.pv()].
#' @param ...     further arguments passed to or from other methods.
#'
#' @details
#' Simply returns the results of the `get_summary_class()` method of the
#' underlying `object`.
#'
#' @returns
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
