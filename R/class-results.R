#' @title
#' Discrete Test Results Class
#'
#' @description
#' This is the class used by \code{DiscreteTests} to output more detailed
#' results, if the \code{simple.output} parameter is set to \code{FALSE}. All
#' data members are private to prevent causing inconsistencies by deliberate or
#' inadvertent changes. However, the results can be read by public methods.
#'
#'
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_character assert_numeric assert_choice assert_class assert_integerish assert_list
#' @export
DiscreteTestResults <- R6Class(
  "DiscreteTestResults",

  ## public ----

  public = list(
    #' @description
    #' Creates a new \code{DiscreteTestResults} object.
    #'
    #' @param test_name        a single character string with the name of the
    #'                         test(s).
    #' @param parameters       a named list containing vectors of the tests
    #'                         parameters. The first element is always a vector
    #'                         of the observed values, all other list items
    #'                         contain the UNIQE parameter combinations of the
    #'                         tests.
    #' @param alternative      a single character string with the testing
    #'                         alternative, e.g. "two.sided".
    #' @param p_values         a numeric vector of the p-values calculated by
    #'                         each discrete testing scenario.
    #' @param support_values   a list of UNIQUE numeric vectors containing all
    #'                         p-values the respective discrete test setting can
    #'                         produce.
    #' @param support_indices  a list of numeric vectors containing the test
    #'                         indices that indicates to which individual
    #'                         testing scenario each unique support belongs.
    #' @param data_name        a single character string with the name of the
    #'                         variable that contains the observed data.
    initialize = function(
      test_name,
      parameters,
      alternative,
      p_values,
      support_values,
      support_indices,
      data_name
    ) {
      # make sure that test name is a single character string
      assert_character(
        x = test_name,
        any.missing = FALSE,
        len = 1,
        null.ok = FALSE
      )

      # make sure that parameters are given in a list of numeric vectors
      assert_list(
        x = parameters,
        types = c("numeric", "character", "atomicvector"),
        any.missing = FALSE,
        null.ok = TRUE
      )

      # make sure observations are numeric
      assert_numeric(
        parameters[[1]],
        any.missing = FALSE,
        min.len = 1,
        null.ok = FALSE
      )

      # overall number of tests, i.e. observations that were tested
      len <- length(parameters[[1]])

      # make sure all list elements (except observations) have the same length
      if(length(unique(sapply(parameters[-1], length))) > 1)
        stop("All vectors of 'parameters' list, except the first, must have the same length")

      # make sure all parameters have correct lengths
      for(i in 1 + seq_along(parameters[-1])){
        if(is.numeric(parameters[[i]])){
          assert_numeric(
            x = parameters[[i]],
            any.missing = FALSE
          )
        }else{
          assert_character(
            x = parameters[[i]],
            any.missing = FALSE,
            len = len
          )
        }
      }

      # make sure that test alternative is a single character string
      assert_character(
        x = alternative,
        any.missing = FALSE,
        min.len = 1,
        max.len = len
      )

      # make sure alternative is one of the 'usual'  ones
      assert_choice(
        x = alternative,
        choices = c("greater", "less", "two.sided", "minlike",
                    "blaker", "central", "absdist"),
        null.ok = FALSE
      )

      # make sure that vector of p-values is numeric with probabilities in [0, 1]
      assert_numeric(
        x = p_values,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        len = len
      )

      # make sure that list of support values is a list
      assert_list(
        x = support_values,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1,
        null.ok = FALSE
      )

      for(i in seq_along(support_values)){
        # make sure each list item contains vectors of probabilities in [0, 1]
        assert_numeric(
          x = support_values[[i]],
          lower = 0,
          upper = 1,
          any.missing = FALSE,
          min.len = 1,
          sorted = TRUE,
          null.ok = FALSE
        )
      }

      # set of p-value indices for checking of correct indices in supports list
      idx_set <- 1L:len

      # make sure that list of support indices is a list
      assert_list(
        x = support_indices,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1,
        null.ok = FALSE
      )

      for(i in seq_along(support_indices)){
       # make sure indices are integerish vectors (and coerce them)
        support_indices[[i]] <- assert_integerish(
          x = support_indices[[i]],
          lower = 1,
          upper = len,
          any.missing = FALSE,
          min.len = 1,
          unique = TRUE,
          coerce = TRUE,
          null.ok = FALSE
        )

        # remove indices of current list item from check set
        idx_set <- setdiff(idx_set, support_indices[[i]])
      }

      # make sure check set is empty
      if(length(idx_set))
        stop("All support set indices must be unique")

      # make sure that data variable name is a single character string
      assert_character(
        x = data_name,
        len = 1,
        null.ok = TRUE
      )

      # assign inputs
      private$test_name         <- test_name
      private$parameters        <- parameters
      private$alternative       <- alternative
      private$p_values          <- p_values
      private$support_values    <- support_values
      private$support_indices   <- support_indices
      private$data_name         <- data_name
    },

    #' @description
    #' Returns the computed p-values.
    get_pvalues = function(){
      return(private$p_values)
    },

    #' @description
    #' Returns the unique test supports.
    get_support_values = function(){
      return(private$support_values)
    },

    #' @description
    #' Returns the indices that indicate to which testing scenario each
    #' unique support belongs.
    get_support_indices = function(){
      return(private$support_indices)
    },

    #' @description
    #' Return the list of the test parameters.
    get_parameters = function(){
      return(private$parameters)
    },

    #' @description
    #' Returns the summary table which is printed by \code{print()}. For better
    #' readability, the packages functions pass non-syntactic parameter names.
    #' As a result, the returned summary table may have non syntactic column
    #' names.
    summary = function(){
      if(is.null(private$summary_table) || !nrow(private$summary_table)){#len <- length(private$p_values)
        idx_scns <- unlist(private$support_indices)
        idx_lens <- sapply(private$support_indices, length)
        par_tab <- sapply(private$parameters[-1], rep, idx_lens)[order(idx_scns), ]
        df <- data.frame(
          private$parameters[[1]],
          par_tab,
          #private$parameters,
          private$p_values
        )
        #names(df) <- c(private$statistics_name, names(private$parameters), "P-Value")
        names(df) <- c(names(private$parameters), "P-Value")
        private$summary_table <- df#as_tibble(df)
      }
      private$summary_table
    },

    #' @description
    #' Returns the computed p-values.
    #' @param ...  further arguments passed to \code{print.data.frame}.
    print = function(...){
      header <- paste0(
        ifelse(!(private$alternative %in% c("greater", "less")), "Two", "One"),
        "-sided ",
        private$test_name
      )

      cat("\n")
      cat("\t", header, "\n")
      cat("\n")
      cat("data:  ", private$data_name, "\n", sep = "")
      if(private$alternative %in% c("greater", "less")){
        meth <- ifelse(private$alternative == "less", "lower", "upper")
        cat("One-sided p-values:  ", meth, " tail\n", sep = "")
      }else{
        meth <- switch(
          private$alternative,
          minlike = "minimum likelihood",
          blaker  = "combined tails",
          absdist = "absolute distance from mean",
          "minimum tail doubling"
        )
        cat("Exact two-sided p-values:  ", meth, "\n", sep = "")
      }
      cat("\n")
      print(self$summary(), ...)
      cat("\n")
      self
    }
  ),

  ## private ----

  private = list(
    # numeric vector of the p-values calculated by each discrete test scenario
    p_values          = numeric(),

    # list of UNIQUE numeric vectors containing all p-values a discrete test
    # setting can produce
    support_values    = list(),

    # list of numeric vectors containing the test indices that indicate to
    # which individual test each unique support belongs
    support_indices   = list(),

    # single character string with the name of the test(s)
    test_name         = character(),

    # named list containing the tests parameters
    parameters        = list(),

    # single character string with the testing alternative, e.g. "two.sided"
    alternative       = character(),

    # single character string with the name of the variable that contains the
    # observed data
    data_name         = character(),

    # data frame that summarized the results of all tests
    summary_table     = data.frame()
  )
)

#' @title Summarizing Discrete Test Results
#'
#' @description
#' \code{summary} method for class \code{"DiscreteTestResults"}.
#'
#' @param object  an object of class \code{"DiscreteTestResults"}, usually
#'                produced by a call to one of the packages test functions, e.g.
#'                \code{\link{binom.test.pv}}.
#' @param ...     further arguments passed to or from other methods.
#'
#' @details
#' Simply returns the results of the \code{get_summary_class()} method of the
#' underlying \code{object}.
#'
#' @examples
#' obj <- binom.test.pv(0:5, 5, 0.5, simple.output = FALSE)
#' summary(obj)
#'
#' @importFrom checkmate assert_class
#' @export
## S3 method for class 'DiscreteTestResults'
summary.DiscreteTestResults <- function(object, ...){
  assert_class(object, c("R6", "DiscreteTestResults"))

  object$summary(...)
}
