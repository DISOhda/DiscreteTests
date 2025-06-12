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
#' # binomial tests
#' obj <- binom.test.pv(0:5, 5, 0.5)
#' # create DiscreteTestResultsSummary object
#' res <- DiscreteTestResultsSummary$new(obj)
#' # print summary
#' print(res)
#' # extract summary table
#' res$get_summary_table()
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_r6
#' @importFrom tibble add_column as_tibble rowid_to_column
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
      assert_r6(test_results, "DiscreteTestResults")

      # get test results
      inputs <- test_results$get_inputs(unique = FALSE)
      pvals  <- test_results$get_pvalues(named = TRUE)

      # compile data as list, ensure desired order
      summary_table <- list(
        inputs$observations,
        c(
          as.list(inputs$parameters),
          as.list(inputs$nullvalues),
          as.list(inputs$computation)
        ),
        pvals
      )

      # get test designations (a.k.a names) if present
      test_names <- names(summary_table[[3]])

      # flatten list
      summary_table <- c(
        summary_table[1],
        summary_table[[2]],
        summary_table[3]
      )

      # if observations are samples, set names (and re-arrange if multi-samples)
      if(is.list(summary_table[[1]]) && !is.data.frame(summary_table[[1]])) {
        # determine if one- or multi-sample test
        len_obs <- length(summary_table[[1]])
        if(len_obs > 1L && is.list(summary_table[[1]][[1]])) {
          # multi-samples
          names_obs <- paste("sample", seq_len(len_obs))
          summary_table <- c(
            summary_table[[1]],
            summary_table[seq_along(summary_table)[-1]]
          )
        } else {
          # single samples
          names_obs <- "sample"
        }
      } else {
        # single observations
        names_obs <- names(summary_table[[1]])
        summary_table <- c(
          as.list(summary_table[[1]]),
          summary_table[seq_along(summary_table)[-1]]
        )
      }

      # create tibble
      summary_table <- as_tibble(summary_table, .name_repair = "minimal")

      # set column headers
      names(summary_table) <- c(
        names_obs,
        names(inputs$parameters),
        names(inputs$nullvalues),
        names(inputs$computation),
        "p-value"
      )

      # add ID column
      summary_table <- if(!is.null(test_names))
        add_column(summary_table, ID = test_names, .before = 1) else
          rowid_to_column(summary_table, "ID")

      # assign inputs
      private$test_results  <- test_results
      private$summary_table <- summary_table
    },

    #' @description
    #' Returns the underlying [DiscreteTestResults] object.
    #'
    #' @return
    #' A [DiscreteTestResults] R6 class object.
    get_test_results = function() {
      return(private$test_results)
    },

    #' @description
    #' Returns the summary table of the underlying [DiscreteTestResults] object.
    #' @return
    #' A data frame.
    get_summary_table = function() {
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
    print = function(...) {
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
#' # binomial tests
#' obj <- binom.test.pv(0:5, 5, 0.5)
#' # print summary
#' summary(obj)
#' # extract summary table
#' smry <- summary(obj)
#' smry$get_summary_table()
#'
#' @importFrom checkmate assert_r6
#' @export
## S3 method for class 'DiscreteTestResults'
summary.DiscreteTestResults <- function(object, ...) {
  assert_r6(object, "DiscreteTestResults")

  DiscreteTestResultsSummary$new(object)
}
