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
#'       # mandatory exactness information, needs to be replicated to length of 'x'
#'       exact = rep(TRUE, m),
#'       # distribution information, needs to be replicated to length of 'x'
#'       distribution = "binomial",
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
#'   # test statistics (not needed, since observation itself is the statistic)
#'   statistics = NULL,
#'   # numerical vector of p-values
#'   p_values = pv,
#'   # list of supports (here: only one support); values must be sorted and unique
#'   pvalue_supports = list(unique(sort(support))),
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
#' @importFrom checkmate qassert
#' @importFrom checkmate assert_character assert_subset
#' @importFrom checkmate assert_integerish assert_numeric
#' @importFrom checkmate assert_logical
#' @importFrom checkmate assert_data_frame assert_list
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
    #'                          hypothesised null values **as data frames or**
    #'                          **lists**; the names of these list fields must
    #'                          be `observations`, `nullvalues` and
    #'                          `parameters`. See details for further
    #'                          information about the requirements for these
    #'                          fields.
    #' @param statistics        data frame containing the tests' statistics;
    #'                          `NULL` is allowed and recommended if the
    #'                          observed values themselves are the statistics.
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
    #'   \item{`$observations`}{data frame or list of vectors that contains the
    #'                          observed data; if the observed data is a matrix,
    #'                          it must be converted to a data frame; must not
    #'                          be `NULL`, only numerical and character values
    #'                          are allowed.}
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
    #'                        **must** include mandatory columns named `exact`,
    #'                        distribution and `alternative`; only numerical,
    #'                        character or logical values are allowed.}
    #' }
    #'
    #' Missing values or `NULL`s are not allowed for any of these fields. All
    #' data frames must have the same number of rows. Their column names are
    #' used by the `print()` method for producing text output, therefore they
    #' should be informative, i.e. short and (if necessary) non-syntactic,
    #' like e.g. `` `number of success` ``.
    #'
    #' The mandatory column `exact` of the data frame `parameters` must be
    #' logical, while `distribution` must contain a character string describing
    #' the distribution under the null hypothesis, e.g. `"normal"`, and is used
    #' by the `print()` method. The column `alternative` must contain one of the
    #' strings `"greater"`, `"less"`, `"two.sided"`, `"minlike"`, `"blaker"`,
    #' `"absdist"` or `"central"`.
    #'
    initialize = function(
      test_name,
      inputs,
      statistics,
      p_values,
      pvalue_supports,
      support_indices,
      data_name
    ) {
      # ensure that test name is a single character string
      qassert(x = test_name, rules = "S1")

      # ensure that inputs are given as data frames or lists
      assert_list(
        x = inputs,
        types = c("data.frame", "list"),
        len = 3L,
        names = "named"
      )

      # ensure that input list elements are properly named
      if(any(!(c("observations", "parameters", "nullvalues") %in% names(inputs))))
        stop("Names of fields in list 'inputs' must be 'observations', 'parameters' and 'nullvalues'")

      # ensure that observations are in a proper format depending on type
      if(is.data.frame(inputs$observations)) {
        # ensure observations are in a data frame of numerical, character or
        #   logical vectors
        assert_data_frame(
          x = inputs$observations,
          types = c("numeric", "character", "logical"),
          any.missing = FALSE,
          min.cols = 1L,
          min.rows = 1L
        )
        # overall number of tests, i.e. observations that were tested
        len <- nrow(inputs$observations)
      } else {
        # ensure observations are in a list containing numerical or character
        #   vectors
        assert_list(
          x = inputs$observations,
          types = c("numeric", "character", "logical", "list"),
          any.missing = FALSE,
          min.len = 1L
        )
        # if observations are a list, they must contain sample tuples
        if(is.list(inputs$observations)) {
          # all list elements must have the same data type
          types <- unique(sapply(inputs$observations, mode))
          if(length(types) > 1L)
            stop("All observations must have the same data type")
          # if type is list: multiple samples; else: single samples
          if(types == "list") {
            len <- unique(sapply(inputs$observations, length))
            if(length(len) > 1L)
              stop("All sample lists must have the same length")
            # for(i in seq_len(len)) {
            #   types <- unique(sapply(inputs$observations[[i]], mode))
            #   if(length(types) > 1L)
            #     stop("All samples must have the same data type")
            #   if(!(types %in% c("numeric", "character", "logical")))
            #     stop("All samples must be numeric, character or logical")
            # }
            types <- character(0)
            for(i in seq_along(inputs$observations))
              types <- unique(c(types, sapply(inputs$observations[[i]], mode)))
            if(length(types) > 1L)
              stop("All samples must have the same data type")
          } else {
            # overall number of tests that were performed
            len <- length(inputs$observations)
          }
          if(!(types %in% c("numeric", "character", "logical")))
            stop("All samples must be numeric, character or logical")
        }
      }

      # ensure that the parameters are in a data frame with as many rows as
      #   there are observations and that it contains only numerical, character
      #   or logical vectors
      assert_data_frame(
        x = inputs$parameters,
        types = c("numeric", "character", "logical"),
        any.missing = TRUE,
        nrows = len
      )

      # ensure that alternatives exist and that they are strings
      if(exists("alternative", inputs$parameters)) {
        assert_character(
          x = inputs$parameters$alternative,
          any.missing = FALSE
        )
      } else stop("'parameters' must have an 'alternative' column")

      # ensure that information about exactness exists and that it is logical
      if(exists("exact", inputs$parameters)) {
        assert_logical(
          x = inputs$parameters$exact,
          any.missing = FALSE
        )
      } else stop("'parameters' must have an 'exact' column")

      # ensure that information about the null distribution exists and that it
      #   is text
      if(exists("distribution", inputs$parameters)) {
        assert_character(
          x = inputs$parameters$distribution,
          any.missing = FALSE
        )
      } else stop("'parameters' must have an 'distribution' column")

      # ensure that all alternatives are from the 'usual'  ones
      assert_subset(
        x = inputs$parameters$alternative,
        choices = c("greater", "less", "two.sided", "minlike",
                    "blaker", "central", "absdist"),
        empty.ok = FALSE
      )

      # ensure that all null (i.e. hypothesised) values are in a numeric data
      #   frame
      assert_data_frame(
        x = inputs$nullvalues,
        types = "numeric",
        any.missing = FALSE,
        nrows = len
      )

      # ensure that the statistics are in a data frame that contains only
      #   numerical vectors and that it has as many rows as there are
      #   observations (NULL is allowed if there are no statistics)
      assert_data_frame(
        x = statistics,
        types = "numeric",
        any.missing = FALSE,
        nrows = len,
        null.ok = TRUE
      )

      # ensure that vector of p-values is numeric with probabilities in [0, 1]
      assert_numeric(
        x = p_values,
        lower = 0L,
        upper = 1L,
        any.missing = FALSE,
        len = len
      )

      # ensure that list of support values is a list
      assert_list(
        x = pvalue_supports,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1L,
        max.len = len
      )

      for(i in seq_along(pvalue_supports)){
        # ensure each list item contains sorted vectors of probabilities in
        #   [0, 1]
        assert_numeric(
          x = pvalue_supports[[i]],
          lower = 0L,
          upper = 1L,
          any.missing = FALSE,
          min.len = 1L,
          sorted = TRUE
        )
      }

      # set of p-value indices for checking of correct indices in supports list
      idx_set <- seq_len(len)

      # ensure that list of support indices is a list
      assert_list(
        x = support_indices,
        types = "numeric",
        any.missing = FALSE,
        min.len = 1L
      )

      for(i in seq_along(support_indices)){
       # ensure indices are integerish vectors (and coerce them)
        support_indices[[i]] <- assert_integerish(
          x = support_indices[[i]],
          lower = 1L,
          upper = len,
          any.missing = FALSE,
          min.len = 1L,
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
      private$statistics      <- statistics
      private$p_values        <- p_values
      private$pvalue_supports <- pvalue_supports
      private$support_indices <- support_indices
      private$data_name       <- data_name

      # add row names in input data (if present) to p-values
      names(private$p_values) <- rownames(inputs$observations)
    },

    #' @description
    #' Returns the computed p-values.
    #'
    #' @param named  single logical value that indicates whether the vector is
    #'               to be returned as a named vector (if names are present)
    #'
    #' @return
    #' A numeric vector of the p-values of all null hypotheses.
    #'
    get_pvalues = function(named = TRUE) {
      if(named)
        return(private$p_values) else
          return(as.numeric(private$p_values))
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
      lst <- private$inputs[c("observations", "nullvalues", "parameters")]
      if(unique) {
        nc <- ncol(lst$nullvalues)
        df <- unique(cbind(lst$nullvalues, lst$parameters))
        lst$nullvalues <- df[ seq_len(nc)]
        lst$parameters <- df[-seq_len(nc)]
      }
      return(lst)
    },

    #' @description
    #' Returns the test statistics.
    #'
    #' @return
    #' A numeric `data.frame` with one column containing the test statistics.
    #'
    get_statistics = function() {
      return(private$statistics)
    },

    #' @description
    #' Returns the *p*-value supports, i.e. all observable p-values under the
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
      if(!is.null(test_idx) && !all(is.na(test_idx)))
        qassert(test_idx, "x+[0,)")
      qassert(limit, c("0", "x1", "X1[0,)"))
      assert_int(x = limit, na.ok = TRUE, lower = 0, null.ok = TRUE)

      if(inputs || pvalues){
        pars <- self$get_inputs(unique = FALSE)
        idx_alt <- which(names(pars$parameters) == "alternative")
        idx_ex <- which(names(pars$parameters) == "exact")
        idx_dist <- which(names(pars$parameters) == "distribution")
        idx <- c(idx_alt, idx_ex, idx_dist)
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

        names <- if(is.data.frame(pars$observations))
          rownames(pars$observations)

        for(i in nums) {
          # number of digits of i
          chars_i <- floor(log10(i)) + 1
          cat("\n")
          cat("Test", i)
          if(!is.null(names) && names[i] != i) {
            # how many characters can be printed to console (minus "Test", " ", ":")
            chars_tst_line <- getOption("width") - 6
            # how many characters can be used by tag (= row name)
            chars_name_max <- chars_tst_line - 3
            # length of current row name
            chars_name <- nchar(names[i])
            # print row name in brackets behind test number
            if(chars_name > chars_name_max) {
              # abbreviate names that are too long and add "..."
              cat(
                " (",
                substr(gsub("[\r\n]", "_", names[i]), 1, chars_name_max - 3),
                "...)",
                sep = ""
              )
            } else cat(" (", gsub("[\r\n]", "_", names[i]), ")", sep = "")
            # length of names (plus " ", "(" and ")") for underscores
            chars_name <- min(chars_name, chars_name_max) + 3
          } else chars_name <- 0
          cat(":\n")
          for(j in 1:(6 + chars_name + chars_i)) cat("-")
          cat("\n")

          if(inputs) {
            # print function with wrapping and indentation for observations,
            #  parameters and hypotheses
            print_wrap <- function(alt, str, indent) {
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

            # width for all indentations/wrappings
            indent_width <- 18

            # print observations
            # number_samples <- if(is.data.frame(pars$observations)) {
            #   ncol(pars$observations)
            # } else {
            #   ifelse(
            #     is.list(pars$observations[[i]]),
            #     length(pars$observations[[i]]),
            #     1
            #   )
            # }
            # str <- character(number_samples)

            if(is.data.frame(pars$observations)) {
              str <- paste(
                paste0(
                  names(pars$observations), " = ", pars$observations[i, ]
                ),
                collapse = ", "
              )
              print_wrap("observations:", str, indent_width)
            } else {
              number_samples <- ifelse(
                is.list(pars$observations[[1]]),
                length(pars$observations),
                1
              )
              for(s in seq_len(number_samples)) {
                # sample_obs <- pars$observations[[s]]
                # if(number_samples > 1) sample_obs <- sample_obs[[s]]
                sample_obs <- if(number_samples > 1)
                  pars$observations[[s]][[i]] else
                    pars$observations[[i]]
                len_obs <- length(sample_obs)
                str <- paste(
                  len_obs,
                  class(sample_obs),
                  ifelse(len_obs > 1, "values", "value"),
                  paste0(
                    "(",
                    paste(
                      format(sample_obs[seq_len(min(6, len_obs))], TRUE),
                      collapse = ", "
                    ),
                    if(len_obs > 7) ", ..." else
                      ifelse(
                        len_obs == 7,
                        paste0(", ", format(sample_obs[7])),
                        ""
                      )
                    ,
                    ")"
                  )
                )
                print_wrap(
                  ifelse(
                    number_samples == 1,
                    "observations:",
                    paste0("sample ", s, ":")
                  ),
                  str,
                  indent_width
                )
              }
            }

            # print statistics (if any)
            if(!is.null(private$statistics)) {
              str <- paste(
                paste0(
                  names(private$statistics), " = ", private$statistics[i, ]
                )
              )
              print_wrap("statistics:", str, indent_width)
            }

            # print computation method
            str <- ifelse(pars$parameters[[idx_ex]][i], "exact", "approximate")
            str <- if(length(idx_dist)) {
              paste0(
                str,
                ", using ",
                pars$parameters[[idx_dist]][i],
                " distribution"
              )
            }
            print_wrap("computation:", str, indent_width)

            # print parameters
            if(ncol(pars$parameters) > 3) {
              idx_par_valid <- which(!is.na(pars$parameters[i, -idx]))
              if(length(idx_par_valid)) {
                str <- paste(
                  paste0(
                    names(pars$parameters)[-idx][idx_par_valid],
                    " = ",
                    format(pars$parameters[i, -idx][idx_par_valid])
                  ),
                  collapse = ", "
                )
                print_wrap("parameters:", str, indent_width)
              }
            }

            if(pars$parameters[[idx_alt]][i] == "greater") {
              alt <- "greater than"
              null <- "less than or equal to"
            } else if(pars$parameters[[idx_alt]][i] == "less") {
              alt <- "less than"
              null <- "greater than or equal to"
            } else {
              alt <- "not equal to"
              null <- "equal to"
            }

            str <- paste("true", names(pars$nullvalues)[1], "is", null,
                         pars$nullvalues[[1]][i])
            print_wrap("null hypothesis:", str, indent_width)

            str <- paste("true", names(pars$nullvalues)[1], "is", alt,
                         pars$nullvalues[[1]][i])
            print_wrap("alternative:", str, indent_width)
          }

          if(pvalues) {
            meth <- switch(
              EXPR    = pars$parameters[[idx_alt]][i],
              minlike = "(minimum likelihood)",
              blaker  = "(combined tails)",
              absdist = "(absolute distance from mean)",
              central = "(minimum tail doubling)",
              ""
            )
            cat("p-value:          ")
            cat(private$p_values[i], meth, "\n")
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
              "use print parameter 'limit' for more results ]\n"
            )
          )
        if(!is.null(test_idx))
          cat(paste("[", length(nums), "out of", n, "results printed ]\n"))
      }

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

    ## @field statistics         data frame containing the tests' statistics
    statistics = data.frame(),

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
