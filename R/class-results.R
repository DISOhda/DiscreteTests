#' @title
#' Discrete Test Results Class
#'
#' @description
#' This is the class used by `DiscreteTests` for returning more detailed p-value
#' results, if the `simple.output` parameter is set to `FALSE`. All data members
#' are private to prevent causing inconsistencies by deliberate or inadvertent
#' changes. However, the results can be read by public methods.
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_character assert_choice assert_class assert_integerish assert_list assert_numeric qassert
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
      qassert(x = test_name, rules = "S1")

      # make sure that inputs are given in a list of numeric vectors
      assert_list(
        x = inputs,
        types = c("numeric", "vector", "matrix", "list", "null"),
        len = 2
      )

      # make sure observations (first list element) are numeric
      qassert(x = inputs[[1]], rules = "N+")

      # overall number of tests, i.e. observations that were tested
      len <- ifelse(
        test = is.matrix(inputs[[1]]) || is.data.frame(inputs[[1]]),
        yes = nrow(inputs[[1]]),
        no = length(inputs[[1]])
      )

      # make sure second element of input list is a list, too
      assert_list(
        x = inputs[[2]],
        types = c("numeric", "character", "atomicvector"),
        any.missing = FALSE,
        null.ok = TRUE
      )

      # make sure all parameters (second list element, which must be a list) have the same length
      if(length(unique(sapply(inputs[[2]], length))) > 1)
        stop("All vectors of second 'inputs' list element must have the same length")

      # make sure all parameter inputs have correct types
      for(i in seq_along(inputs[[2]])){
        qassert(x = inputs[[2]][[i]], rules = c("N+", "S+"))
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
        # make sure each list item contains sorted vectors of probabilities in [0, 1]
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
      qassert(x = data_name, c("S+", "0"))

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
    #' @returns
    #' A numeric vector of the p-values of each null hypothesis.
    get_pvalues = function(){
      return(private$p_values)
    },

    #' @description
    #' Return the list of the test inputs. It can be chosen, if only unique
    #' parameter sets are needed.
    #' @param unique   integer value that indicates whether only unique
    #'                 parameter sets are to be returned. If `unique = FALSE`
    #'                 (the default), the returned supports may be duplicated.
    #' @returns
    #' A list of two elements. The first one contains the observations for each
    #' tested null hypothesis. The second is another list with the parameter
    #' sets.  If `unique = TRUE`, only unique parameter sets are returned.
    get_inputs = function(unique = FALSE){
      if(!unique){
        idx_scns <- unlist(private$scenario_indices)
        idx_lens <- sapply(private$scenario_indices, length)
        lst <- private$inputs
        lst[[2]] <- lapply(
          lst[[2]],
          function(v) rep(v, idx_lens)[order(idx_scns)]
        )

        return(lst)
      }else return(private$inputs)
    },

    #' @description
    #' Returns the testing scenario supports. It can be chosen, if only unique
    #' supports are needed.
    #' @param unique   integer value that indicates whether only unique supports
    #'                 are to be returned. If `unique = FALSE` (the default),
    #'                 the returned supports may be duplicated.
    #' @returns
    #' A list of numeric vectors. Each one contains all observable p-values of
    #' the respective null hypothesis. If `unique = TRUE`, only unique supports
    #' are returned.
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
    #' @returns
    #' A list of numeric vectors. Each one contains the indices of the null
    #' hypotheses to which the respective support and/or parameter set belongs.
    get_scenario_indices = function(){
      return(private$scenario_indices)
    },

    #' @description
    #' Prints the computed p-values.
    #' @param pvalues    a single logical value that indicates if the resulting
    #'                   p-values are to be printed.
    #' @param inputs     a single logical value that indicates if the inputs
    #'                   values (i.e. observations and parameters) are to be
    #'                   printed.
    #' @param supports   a single logical value that indicates if the p-value
    #'                   supports are to be printed.
    #' @param ...        further arguments passed to `print.default`.
    #' @returns
    #' Prints a summary of the tested null hypotheses. The object itself is
    #' invisibly returned.
    print = function(pvalues = TRUE, inputs = TRUE, supports = TRUE, ...){
      qassert(pvalues, "B1")
      qassert(inputs, "B1")
      qassert(supports, "B1")

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
      if(pvalues){
        cat("\n")
        cat("Resulting p-values:\n")
        print(private$p_values, ...)
        cat("\n")
      }
      if(inputs){
        cat("\n")
        cat("Inputs:\n")
        print(self$get_inputs(unique = FALSE), ...)
      }
      if(supports){
        if(!inputs) cat("\n")
        cat("Supports:\n")
        print(self$get_scenario_supports(unique = FALSE), ...)
      }

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
    summary_table = data.frame(),

    # version of class definition
    class_version = "0.1.0"
  )
)
