#' @title
#' Discrete Test Results Class
#'
#' @description
#' This is the class used by the statistical test functions of this package for
#' returning not only p-values, but also the supports of their distributions and
#' the parameters of the respective tests. Objects of this class are obtained by
#' setting the `simple.output` parameter of a test function to `FALSE` (the
#' default). All data members of this class are private to avoid inconsistencies
#' by deliberate or inadvertent changes by the user. However, the results can be
#' read by public methods.
#'
#' @examples
#' ## one-sided binomial test
#' #  parameters
#' x <- 2:4
#' n <- 5
#' p <- 0.4
#' m <- length(x)
#' #  support (same for all three tests) and p-values
#' support <- sapply(0:n, function(k) binom.test(k, n, p, "greater")$p.value)
#' pv <- support[x + 1]
#' #  DiscreteTestResults object
#' res <- DiscreteTestResults$new(
#'   # string with name of the test
#'   test_name = "Exact binomial test",
#'   # list of data frames
#'   inputs = list(
#'     observations = data.frame(
#'       `number of successes` = x,
#'       # no name check of column header to have a speaking name for 'print'
#'       check.names = FALSE
#'     ),
#'     parameters = data.frame(
#'       # parameter 'n', needs to be replicated to length of 'x'
#'       `number of trials` = rep(n, m),
#'       # mandatory parameter 'alternative', needs to be replicated to length of 'x'
#'       alternative = rep("greater", m),
#'       # no name check of column header to have a speaking name for 'print'
#'       check.names = FALSE
#'     ),
#'     nullvalues = data.frame(
#'       # here: only one null value, 'p'; needs to be replicated to length of 'x'
#'       `probability of success` = rep(p, m),
#'       # no name check of column header to have a speaking name for 'print'
#'       check.names = FALSE
#'     )
#'   ),
#'   # numerical vector of p-values
#'   p_values = pv,
#'   # list of supports (here: only one support); values must be sorted and unique
#'   pvalue_supports = list(sort(unique(support))),
#'   # list of indices that indicate which p-value/hypothesis each support belongs to
#'   support_indices = list(1:m),
#'   # name of input data variables
#'   data_name = "x, n and p"
#' )
#'
#' #  print results without supports
#' print(res)
#' #  print results with supports
#' print(res, supports = TRUE)
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_character assert_data_frame assert_integerish assert_list assert_numeric assert_subset qassert
#' @export
DiscreteTestResults <- R6Class(
  "DiscreteTestResults",

  ## public ----

  public = list(
    #' @description
    #' Creates a new `DiscreteTestResults` object.
    #'
    #' @param test_name         single character string with the name of the
    #'                          test(s).
    #' @param inputs            named list of **exactly three** elements
    #'                          containing the observations, test parameters and
    #'                          hypothesised null values **as data frames**;
    #'                          names of these list fields must be
    #'                          `observations`, `nullvalues` and `parameters`.
    #'                          See details for further information about the
    #'                          requirements for these fields.
    #' @param p_values          numeric vector of the p-values calculated by
    #'                          each hypothesis test.
    #' @param pvalue_supports   list of **unique** numeric vectors containing
    #'                          all p-values that are observable under the
    #'                          respective hypothesis; each value of `p_values`
    #'                          must occur in its respective p-value support.
    #' @param support_indices   list of numeric vectors containing the test
    #'                          indices that indicates to which individual
    #'                          testing scenario each unique parameter set and
    #'                          each unique support belongs.
    #' @param data_name         single character string with the name of the
    #'                          variable that contains the observed data.
    #'
    #' @details
    #' The fields of the `inputs` have the following requirements:
    #' \describe{
    #'   \item{`$observations`}{data frame that contains the observed data; if
    #'                          the observed data is a matrix, it must be
    #'                          converted to a data frame; must not be `NULL`,
    #'                          only numerical and character values are
    #'                          allowed.}
    #'   \item{`$nullvalues`}{data frame that contains the hypothesised values
    #'                        of the tests, e.g. the rate parameters for Poisson
    #'                        tests; must not be `NULL`, only numerical values
    #'                        are allowed.}
    #'   \item{`$parameters`}{data frame that holds the parameter combinations
    #'                        of the null distribution of each test (e.g.
    #'                        numbers of Bernoulli trials for binomial tests, or
    #'                        `m`, `n` and `k` for the hypergeometric
    #'                        distribution used by Fisher's Exact Test, which
    #'                        have to be  derived from the observations first);
    #'                        **must** include a mandatory column named
    #'                        `alternative`; only numerical and character values
    #'                        are allowed.}
    #' }
    #'
    #' Missing values or `NULL`s are not allowed for any of these fields. All
    #' data frames must have the same number of rows. Their column names are
    #' used by the `print` method for producing text output, therefore they
    #' should be informative, i.e. short and (if necessary) non-syntactic,
    #' like e.g. `` `number of success` ``.
    #'
    initialize = function(
      test_name,
      inputs,
      p_values,
      pvalue_supports,
      support_indices,
      data_name
    ) {
      # ensure that test name is a single character string
      qassert(x = test_name, rules = "S1")

      # ensure that inputs are given as lists, data frame or numeric/string vectors
      assert_list(
        x = inputs,
        types = c("data.frame"),
        any.missing = TRUE,
        len = 3,
        names = "named"
      )

      # ensure input list elements are named
      # 'observations', 'parameters', 'null values' and 'alternatives'
      if(any(!(c("observations", "parameters", "nullvalues") %in% names(inputs))))
        stop("Names of list 'inputs' must be 'observations', 'parameters' and 'nullvalues'")

      # ensure observations are a data frame with as many rows as p-values
      assert_data_frame(
        x = inputs$observations,
        types = c("numeric", "character"),
        any.missing = FALSE
      )

      # overall number of tests, i.e. observations that were tested
      len <- nrow(inputs$observations)

      # ensure that the parameters are in a data.frame with at least one row,
      #   containing only numbers or strings
      assert_data_frame(
        x = inputs$parameters,
        types = c("numeric", "character"),
        any.missing = FALSE,
        nrows = len
      )

      # ensure that alternatives exist and that they are strings
      if(exists("alternative", inputs$parameters)) {
        assert_character(
          x = inputs$parameters$alternative,
          any.missing = FALSE
        )
      } else stop("'parameters' must have an 'alternative' column")

      # ensure that all alternatives are from the 'usual'  ones
      assert_subset(
        x = inputs$parameters$alternative,
        choices = c("greater", "less", "two.sided", "minlike",
                    "blaker", "central", "absdist"),
        empty.ok = FALSE
      )

      # ensure that all null (i.e. hypothesised) values are in a numeric data frame
      assert_data_frame(
        x = inputs$nullvalues,
        types = "numeric",
        any.missing = FALSE,
        nrows = len
      )

      # ensure that vector of p-values is numeric with probabilities in [0, 1]
      assert_numeric(
        x = p_values,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        len = len
      )

      # ensure that list of support values is a list
      assert_list(
        x = pvalue_supports,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1,
        max.len = len
      )

      for(i in seq_along(pvalue_supports)){
        # ensure each list item contains sorted vectors of probabilities in [0, 1]
        assert_numeric(
          x = pvalue_supports[[i]],
          lower = 0,
          upper = 1,
          any.missing = FALSE,
          min.len = 1,
          sorted = TRUE
        )
      }

      # set of p-value indices for checking of correct indices in supports list
      idx_set <- 1L:len

      # ensure that list of support indices is a list
      assert_list(
        x = support_indices,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1
      )

      for(i in seq_along(support_indices)){
       # ensure indices are integerish vectors (and coerce them)
        support_indices[[i]] <- assert_integerish(
          x = support_indices[[i]],
          lower = 1,
          upper = len,
          any.missing = FALSE,
          min.len = 1,
          unique = TRUE,
          coerce = TRUE
        )

        # remove indices of current list item from check set
        idx_set <- setdiff(idx_set, support_indices[[i]])

        # ensure that each p-value is taken from its respective distribution
        if(!all(p_values[support_indices[[i]]] %in% c(0, pvalue_supports[[i]])))
          stop("All observed p-values must occur in their respective distribution")
      }


      # ensure check set is empty
      if(length(idx_set))
        stop("All support set indices must be unique")

      # ensure that data variable name is a single character string
      qassert(x = data_name, c("S+", "0"))

      # assign inputs
      private$test_name       <- test_name
      private$inputs          <- inputs
      private$p_values        <- p_values
      private$pvalue_supports <- pvalue_supports
      private$support_indices <- support_indices
      private$data_name       <- data_name
    },

    #' @description
    #' Returns the computed p-values.
    #'
    #' @return
    #' A numeric vector of the p-values of all null hypotheses.
    #'
    get_pvalues = function(){
      return(private$p_values)
    },

    #' @description
    #' Return the list of the test inputs.
    #'
    #' @param unique   single logical value that indicates whether only unique
    #'                 combinations of parameter sets and null values are to be
    #'                 returned. If `unique = FALSE` (the default), the returned
    #'                 data frames may contain duplicate sets.
    #'
    #' @return
    #' A list of three elements. The first one contains a data frame with the
    #' observations for each tested null hypothesis, while the second is another
    #' data frame with the hypothesised null values (e.g. `p` for binomial
    #' tests). The third list field holds the parameter sets (e.g. `n` in case
    #' of a binomial test). If `unique = TRUE`, only unique combinations of
    #' parameter sets and null values are returned, but observations remain
    #' unchanged.
    #'
    get_inputs = function(unique = FALSE) {
      if(unique) {
        lst <- private$inputs[c("observations", "nullvalues", "parameters")]
        nc <- ncol(lst$nullvalues)
        df <- unique(cbind(lst$nullvalues, lst$parameters))
        lst$nullvalues <- df[ seq_len(nc)]
        lst$parameters <- df[-seq_len(nc)]
        return(lst)
      } else return(private$inputs)
    },

    #' @description
    #' Returns the p-value supports, i.e. all observable p-values under the
    #' respective null hypothesis of each test.
    #'
    #' @param unique   single logical value that indicates whether only unique
    #'                 p-value supports are to be returned. If `unique = FALSE`
    #'                 (the default), the returned supports may be duplicated.
    #'
    #' @return
    #' A list of numeric vectors containing the supports of the p-value null
    #' distributions.
    #'
    get_pvalue_supports = function(unique = FALSE) {
      if(!unique) {
        idx_scns <- unlist(private$support_indices)
        idx_lens <- sapply(private$support_indices, length)
        return(rep(private$pvalue_supports, idx_lens)[order(idx_scns)])
      } else return(private$pvalue_supports)
    },

    #' @description
    #' Returns the indices that indicate to which testing scenario each
    #' unique support belongs.
    #'
    #' @return
    #' A list of numeric vectors. Each one contains the indices of the null
    #' hypotheses to which the respective support and/or unique parameter set
    #' belongs.
    #'
    get_support_indices = function(){
      return(private$support_indices)
    },

    #' @description
    #' Prints the computed p-values.
    #'
    #' @param inputs     single logical value that indicates if the inputs
    #'                   values (i.e. observations and parameters) are to be
    #'                   printed; defaults to `TRUE`.
    #' @param pvalues    single logical value that indicates if the resulting
    #'                   p-values are to be printed; defaults to `TRUE`.
    #' @param supports   single logical value that indicates if the p-value
    #'                   supports are to be printed; defaults to `FALSE`.
    #' @param test_idx   integer vector giving the indices of the tests whose
    #'                   results are to be printed; if `NULL` (the default),
    #'                   results of every test up to the index specified by
    #'                   `limit` (see below) are printed
    #' @param limit      single integer that indicates the maximum number of
    #'                   test results to be printed; if `limit = 0`, results of
    #'                   every test are printed; ignored if `test_idx` is not
    #'                   set to `NULL`
    #' @param ...        further arguments passed to `print.default`.
    #'
    #' @return
    #' Prints a summary of the tested null hypotheses. The object itself is
    #' invisibly returned.
    #'
    #' @importFrom checkmate assert_int qassert
    print = function(
      inputs = TRUE,
      pvalues = TRUE,
      supports = FALSE,
      test_idx = NULL,
      limit = 10,
      ...
    ) {
      qassert(inputs, "B1")
      qassert(pvalues, "B1")
      qassert(supports, "B1")
      if(!is.null(test_idx) && !is.na(test_idx))
        qassert(test_idx, "x+[0,)")
      qassert(limit, c("0", "x1", "X1[0,)"))
      assert_int(x = limit, na.ok = TRUE, lower = 0, null.ok = TRUE)

      if(inputs || pvalues){
        pars <- self$get_inputs(unique = FALSE)
        idx <- which(names(pars$parameters) == "alternative")
      }
      if(supports)
        supp <- self$get_pvalue_supports(unique = FALSE)

      cat("\n")
      cat(strwrap(private$test_name, prefix = "\t"), "\n")
      cat("\n")
      cat("data: ", private$data_name, "\n")
      cat("number of tests: ", length(private$p_values), "\n")

      if(any(inputs, pvalues, supports)) {
        n <- length(private$p_values)
        if(is.null(test_idx)) {
          limit <- ifelse(is.null(limit) || is.na(limit), n, limit)
          nums <- seq_len(ifelse(limit, min(limit, n), n))
        } else nums <- unique(pmin(na.omit(test_idx), n))

        for(i in nums) {
          cat("\n")
          cat("Test ", i, ":\n", sep = "")
          for(j in 1:(7 + floor(log10(i)))) cat("-")
          cat("\n")

          if(inputs) {
            # print function with wrapping and indentation for observations,
            #  parameters and hypotheses
            print_wrap <- function(alt, str, indent){
              cat(
                paste0(alt, paste(rep(" ", indent - nchar(alt)), collapse = "")),
                strwrap(
                  x = str,
                  width = getOption("width") - indent,
                  exdent = indent,
                  prefix = "\n",
                  initial = ""
                ),
                "\n",
                sep = ""
              )
            }

            # print observations
            str <- paste(
              paste0(
                names(pars$observations), " = ", pars$observations[i, ]
              ),
              collapse = ", "
            )
            print_wrap("observations:", str, 18)

            # print parameters
            if(ncol(pars$parameters) > 1){
              str <- paste(
                paste0(
                  names(pars$parameters)[-idx], " = ", pars$parameters[i, -idx]
                ),
                collapse = ", "
              )
              print_wrap("parameters:", str, 18)
            }

            if(pars$parameters[[idx]][i] == "greater") {
              alt <- "greater than"
              null <- "less than or equal to"
            } else if(pars$parameters[[idx]][i] == "less") {
              alt <- "less than"
              null <- "greater than or equal to"
            } else {
              alt <- "not equal to"
              null <- "equal to"
            }

            str <- paste("true", names(pars$nullvalues)[1], "is", null,
                         pars$nullvalues[[1]][i])
            print_wrap("null hypothesis:", str, 18)

            str <- paste("true", names(pars$nullvalues)[1], "is", alt,
                         pars$nullvalues[[1]][i])
            print_wrap("alternative:", str, 18)
          }

          if(pvalues) {
            if(pars$parameters[[idx]][i] %in% c("greater", "less")) {
              meth <- switch(
                pars$parameters[[idx]][i],
                greater = "upper tail",
                less = "lower tail"
              )
            } else {
              meth <- switch(
                pars$parameters[[idx]][i],
                minlike = "minimum likelihood",
                blaker  = "combined tails",
                absdist = "absolute distance from mean",
                "minimum tail doubling"
              )
            }
            cat("p-value:          ")
            cat(private$p_values[i], " (", meth, ")\n", sep = "")
          }

          if(supports) {
            cat("support:\n")
            print(supp[[i]], ...)
          }
        }
        cat("\n")
        if(is.null(test_idx) && limit > 0 && limit < n)
          cat(
            paste(
              "[ print limit reached --", n - limit, "results omitted --",
              "use print parameter 'limit' for more results ]"
            )
          )
        if(!is.null(test_idx))
          cat(paste("[", length(nums), "out of", n, "results printed ]"))
      }
      cat("\n")

      self
    }
  ),

  ## private ----

  private = list(
    ## @field test_name          single character string with the name of the
    ##                           test(s), e.g. "Fisher's Exact Test"
    test_name = character(),

    ## @field inputs             named list containing the observations and the
    ##                           tests' parameters
    inputs = list(),

    ## @field p_values           numeric vector of the p-values calculated for
    ##                           each discrete test setting
    p_values = numeric(),

    ## @field pvalue_supports    list of UNIQUE numeric vectors containing all
    ##                           p-values that can be observed in the respective
    ##                           discrete test setting
    pvalue_supports = list(),

    ## @field support_indices    list of numeric vectors containing the test
    ##                           indices that indicate to which individual test
    ##                           setting each unique support belongs
    support_indices = list(),

    ## @field data_name          single character string with the name of the
    ##                           variable that contains the observed data
    data_name = character()
  )
)
