#' @title
#' Discrete Test Results Class
#'
#' @description
#' This is the class used by `DiscreteTests` to output more detailed
#' results, if the `simple.output` parameter is set to `FALSE`. All
#' data members are private to prevent causing inconsistencies by deliberate or
#' inadvertent changes. However, the results can be read by public methods.
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_character assert_numeric assert_choice assert_class assert_integerish assert_list
#' @export
DiscreteTestResults <- R6Class(
  "DiscreteTestResults",

  ## public ----

  public = list(
    #' @description
    #' Creates a new `DiscreteTestResults` object.
    #'
    #' @param test_name           a single character string with the name of the
    #'                            test(s).
    #' @param inputs              a named list of **exactly two** elements
    #'                            containing the tests parameters. The first
    #'                            element holds the observed data, e.g. a vector
    #'                            or a matrix. The second one is a named list of
    #'                            the **unique** parameter combinations of the
    #'                            tests, given by vectors.
    #' @param alternative         a single character string with the testing
    #'                            alternative, e.g. "two.sided".
    #' @param p_values            a numeric vector of the p-values calculated by
    #'                            each discrete testing scenario.
    #' @param scenario_supports   a list of **unique** numeric vectors
    #'                            containing all p-values the respective
    #'                            discrete test setting can produce.
    #' @param scenario_indices    a list of numeric vectors containing the test
    #'                            indices that indicates to which individual
    #'                            testing scenario each unique parameter set and
    #'                            each unique support belongs.
    #' @param data_name           a single character string with the name of the
    #'                            variable that contains the observed data.
    initialize = function(
      test_name,
      inputs,
      alternative,
      p_values,
      scenario_supports,
      scenario_indices,
      data_name
    ) {
      # make sure that test name is a single character string
      assert_character(
        x = test_name,
        any.missing = FALSE,
        len = 1
      )

      # make sure that inputs are given in a list of numeric vectors
      assert_list(
        x = inputs,
        types = c("numeric", "vector", "matrix", "list", "null"),
        len = 2
      )

      # make sure observations (first list element) are numeric
      assert_numeric(
        inputs[[1]],
        any.missing = FALSE,
        min.len = 1
      )

      # overall number of tests, i.e. observations that were tested
      len <- ifelse(is.matrix(inputs[[1]]) || is.data.frame(inputs[[1]]), nrow(inputs[[1]]), length(inputs[[1]]))

      # make sure second element of input list is a list, too
      assert_list(
        x = inputs[[2]],
        types = c("numeric", "character", "atomicvector"),
        any.missing = FALSE,
        null.ok = TRUE
      )

      # make sure all parameters (second list element, which is a list) have the same length
      if(length(unique(sapply(inputs[[2]], length))) > 1)
        stop("All vectors of 'inputs' list, except the first, must have the same length")

      # make sure all parameter inputs have correct lengths
      for(i in seq_along(inputs[[2]])){
        if(is.numeric(inputs[[2]][[i]])){
          assert_numeric(
            x = inputs[[2]][[i]],
            any.missing = FALSE
          )
        }else{
          assert_character(
            x = inputs[[2]][[i]],
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
                    "blaker", "central", "absdist")
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
        x = scenario_supports,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1
      )

      for(i in seq_along(scenario_supports)){
        # make sure each list item contains vectors of probabilities in [0, 1]
        assert_numeric(
          x = scenario_supports[[i]],
          lower = 0,
          upper = 1,
          any.missing = FALSE,
          min.len = 1,
          sorted = TRUE
        )
      }

      # set of p-value indices for checking of correct indices in supports list
      idx_set <- 1L:len

      # make sure that list of support indices is a list
      assert_list(
        x = scenario_indices,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1
      )

      for(i in seq_along(scenario_indices)){
       # make sure indices are integerish vectors (and coerce them)
        scenario_indices[[i]] <- assert_integerish(
          x = scenario_indices[[i]],
          lower = 1,
          upper = len,
          any.missing = FALSE,
          min.len = 1,
          unique = TRUE,
          coerce = TRUE
        )

        # remove indices of current list item from check set
        idx_set <- setdiff(idx_set, scenario_indices[[i]])
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
      private$inputs            <- inputs
      private$alternative       <- alternative
      private$p_values          <- p_values
      private$scenario_supports <- scenario_supports
      private$scenario_indices  <- scenario_indices
      private$data_name         <- data_name
    },

    #' @description
    #' Returns the computed p-values.
    get_pvalues = function(){
      return(private$p_values)
    },

    #' @description
    #' Returns the testing scenario supports. It can be chosen, if only unique
    #' supports are needed.
    #' @param unique   integer value that indicates whether only unique supports
    #'                 are to be returned. If `unique = FALSE` (the default),
    #'                 the returned supports may be duplicated.
    get_scenario_supports = function(unique = FALSE){
      if(!unique){
        idx_scns <- unlist(private$scenario_indices)
        idx_lens <- sapply(private$scenario_indices, length)
        return(rep(private$scenario_supports, idx_lens)[order(idx_scns)])
      } else return(private$scenario_supports)
    },

    #' @description
    #' Returns the indices that indicate to which testing scenario each
    #' unique support belongs.
    get_scenario_indices = function(){
      return(private$scenario_indices)
    },

    #' @description
    #' Return the list of the test inputs
    get_inputs = function(){
      return(private$inputs)
    },

    #' @description
    #' Returns the summary table which is printed by `print()`. For better
    #' readability, the packages functions pass non-syntactic parameter names.
    #' As a result, the returned summary table may have non syntactic column
    #' names.
    summary = function(){
      if(is.null(private$summary_table) || !nrow(private$summary_table)){#len <- length(private$p_values)
        idx_scns <- unlist(private$scenario_indices)
        idx_lens <- sapply(private$scenario_indices, length)
        par_tab <- if(!is.null(private$inputs[[2]]))
          sapply(private$inputs[[2]], rep, idx_lens)[order(idx_scns), ]

        df <- list(
          private$inputs[[1]],
          par_tab,
          private$p_values
        )
        df <- as.data.frame(df[!sapply(df, is.null)])
        if(is.matrix(private$inputs[[1]]) || is.data.frame(private$inputs[[1]]))
          name_obs <- colnames(private$inputs[[1]]) else
            name_obs <- names(private$inputs[1])
        names(df) <- c(name_obs, names(private$inputs[[2]]), "p-value")

        private$summary_table <- df
      }
      private$summary_table
    },

    #' @description
    #' Returns the computed p-values.
    #' @param ...  further arguments passed to `print.data.frame`.
    print = function(...){
      cat("\n")
      cat(strwrap(private$test_name, prefix = "\t"), "\n")
      cat("\n")
      cat("data:  ", private$data_name, "\n", sep = "")
      if(private$alternative %in% c("greater", "less")){
        meth <- switch(
          private$alternative,
          greater = "upper tail",
          less = "lower tail"
        )
        side <- "one"
      }else{
        meth <- switch(
          private$alternative,
          minlike = "minimum likelihood",
          blaker  = "combined tails",
          absdist = "absolute distance from mean",
          "minimum tail doubling"
        )
        side <- "two"
      }
      cat(side, "-sided p-values:  ", meth, "\n", sep = "")
      cat("\n")
      print(self$summary(), ...)
      cat("\n")
      self
    }
  ),

  ## private ----

  private = list(
    # numeric vector of the p-values calculated by each discrete test scenario
    p_values = numeric(),

    # list of UNIQUE numeric vectors containing all p-values a discrete test
    # setting can produce
    scenario_supports = list(),

    # list of numeric vectors containing the test indices that indicate to
    # which individual test each unique support belongs
    scenario_indices = list(),

    # single character string with the name of the test(s)
    test_name = character(),

    # named list containing the tests parameters
    inputs = list(),

    # single character string with the testing alternative, e.g. "two.sided"
    alternative = character(),

    # single character string with the name of the variable that contains the
    # observed data
    data_name = character(),

    # data frame that summarized the results of all tests
    summary_table = data.frame()
  )
)

#' @title Summarizing Discrete Test Results
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
#' @examples
#' obj <- binom.test.pv(0:5, 5, 0.5)
#' summary(obj)
#'
#' @importFrom checkmate assert_class
#' @export
## S3 method for class 'DiscreteTestResults'
summary.DiscreteTestResults <- function(object, ...){
  assert_class(object, c("R6", "DiscreteTestResults"))

  object$summary(...)
}
