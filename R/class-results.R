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
#' @importFrom checkmate assert_int assert_integerish assert_numeric
#' @importFrom checkmate assert_logical
#' @importFrom checkmate assert_data_frame assert_list
#' @importFrom cli cli_end cli_li cli_ul
#' @importFrom cli ansi_collapse cli_h1 cli_h2 cli_h3 cli_text cli_verbatim qty
#' @importFrom cli col_blue col_green col_red
#' @importFrom cli cli_abort
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
        types = c("data.frame", "list", "null"),
        len = 4L,
        names = "named"
      )

      # ensure that input list elements are properly named
      if(!all(
        c("observations", "parameters", "nullvalues", "computation") %in%
          names(inputs)
      ))
        cli_abort(paste(
          "Names of fields in list 'inputs' must be 'observations',",
          "'parameters', 'nullvalues' and 'computation'."
        ))

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
      } else
        if(is.list(inputs$observations)) {
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
            type <- unique(sapply(inputs$observations, mode))
            if(length(type) > 1L)
              cli_abort("All observations must have the same data type")
            # if type is list: multiple samples; else: single samples
            if(type == "list") {
              len <- unique(sapply(inputs$observations, length))
              if(length(len) > 1L)
                cli_abort("All sample lists must have the same length")

              type <- character(0)
              for(i in seq_along(inputs$observations))
                type <- unique(c(type, sapply(inputs$observations[[i]], mode)))
              if(length(type) > 1L)
                cli_abort("All samples must have the same data type")
            } else {
              # overall number of tests that were performed
              len <- length(inputs$observations)
            }
            if(!(type %in% c("numeric", "character", "logical")))
              cli_abort("All samples must be numeric, character or logical")
          }
        } else cli_abort("'observations' must be a list or a data frame")

      # ensure that the parameters are in a data frame with as many rows as
      #   there are observations and that it contains only numerical, character
      #   or logical vectors
      assert_data_frame(
        x = inputs$parameters,
        types = c("numeric", "character", "logical"),
        any.missing = TRUE,
        nrows = len,
        null.ok = TRUE
      )

      # ensure that all null (i.e. hypothesised) values are in a numeric data
      #   frame
      assert_data_frame(
        x = inputs$nullvalues,
        types = "numeric",
        any.missing = FALSE,
        nrows = len,
      )

      # ensure that the computation information is given in a data frame with as
      #   many rows as there are observations and that it contains only
      #   numerical, character or logical vectors
      assert_data_frame(
        x = inputs$computation,
        types = c("numeric", "character", "logical"),
        any.missing = TRUE,
        nrows = len
      )

      # ensure that alternatives exist and that they are strings
      if(exists("alternative", inputs$computation)) {
        assert_character(
          x = inputs$computation$alternative,
          any.missing = FALSE
        )
      } else cli_abort("'computation' must have an 'alternative' column")

      # ensure that all alternatives are from the 'usual' ones
      assert_subset(
        x = inputs$computation$alternative,
        choices = c("greater", "less", "two.sided", "minlike",
                    "blaker", "central", "absdist"),
        empty.ok = FALSE
      )

      # ensure that information about exactness exists and that it is logical
      if(exists("exact", inputs$computation)) {
        assert_logical(
          x = inputs$computation$exact,
          any.missing = FALSE
        )
      } else cli_abort("'computation' must have an 'exact' column")

      # ensure that information about the null distribution exists and that it
      #   is text
      if(exists("distribution", inputs$computation)) {
        assert_character(
          x = inputs$computation$distribution,
          any.missing = FALSE
        )
      } else cli_abort("'computation' must have a 'distribution' column")

      # # ensure that information about potential continuity correction exists and
      # #   that it is logical (NA allowed for exact computation)
      # if(exists("correct", inputs$computation)) {
      #   assert_logical(
      #     x = inputs$computation$correct,
      #     any.missing = TRUE
      #   )
      # } else cli_abort("'computation' must have a 'correct' column")

      # ensure that the statistics are in a data frame that contains only
      #   numerical vectors and that it has as many rows as there are
      #   observations (NULL is allowed if there are no statistics)
      assert_data_frame(
        x = statistics,
        types = "numeric",
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
          cli_abort(paste(
            "All observed p-values must occur in their respective distribution"
          ))
      }


      # ensure check set is empty
      if(length(idx_set))
        cli_abort("All support set indices must be unique")

      # ensure that data variable name is a character vector or NULL
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
    #' A list of four elements. The first one contains a data frame with the
    #' observations for each tested null hypothesis, while the second is another
    #' data frame with additional parameters (if any, e.g. `n` in case of a
    #' binomial test). The third list field holds the hypothesised null values
    #' (e.g. `p` for binomial tests). The last list element contains the
    #' computational parameters, e.g. test `alternative`s, the used
    #' `distribution` etc. If `unique = TRUE`, only unique combinations of
    #' parameters, null values and computation specifics are returned, but
    #' observations remain unchanged (i.e. not unique).
    #'
    get_inputs = function(unique = FALSE) {
      lst <- private$inputs[c(
        "observations", "parameters", "nullvalues", "computation"
      )]
      if(unique) {
        np  <- ncol(lst$parameters)
        nn  <- ncol(lst$nullvalues)
        nc  <- ncol(lst$computation)
        id1 <- seq_len(np)
        id2 <- seq_len(np + nn)
        id3 <- seq_len(np + nn + nc)
        df  <- unique(cbind(lst$parameters, lst$nullvalues, lst$computation))
        lst$parameters  <- df[id1]
        lst$nullvalues  <- df[setdiff(id2, id1)]
        lst$computation <- df[setdiff(id3, id2)]
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
    #' @param inputs            single logical value that indicates if the
    #'                          input values (i.e. observations, statistics and
    #'                          parameters) are to be printed; defaults to
    #'                          `TRUE`.
    #' @param pvalue_details    single logical value that indicates if details
    #'                          about the p-value computation are to be printed;
    #'                          defaults to `TRUE`.
    #' @param supports          single logical value that indicates if the
    #'                          p-value supports are to be printed; defaults to
    #'                          `FALSE`.
    #' @param test_idx          integer vector giving the indices of the tests
    #'                          whose results are to be printed; if `NULL` (the
    #'                          default), results of every test up to the index
    #'                          specified by `limit` (see below) are printed.
    #' @param limit             single integer that indicates the maximum number
    #'                          of test results to be printed; if `limit = 0`,
    #'                          results of every test are printed; ignored if
    #'                          `test_idx` is not set to `NULL`
    #' @param ...               further arguments passed to
    #'                          [`print.default()`][base::print.default()].
    #'
    #' @return
    #' Prints a summary of the tested null hypotheses. The object itself is
    #' invisibly returned.
    #'
    print = function(
      inputs = TRUE,
      pvalue_details = TRUE,
      supports = FALSE,
      test_idx = NULL,
      limit = 10,
      ...
    ) {
      # check arguments for plausibility
      qassert(inputs, "B1")
      qassert(pvalue_details, "B1")
      qassert(supports, "B1")
      if(!is.null(test_idx) && !all(is.na(test_idx)))
        qassert(test_idx, "x+[0,)")
      qassert(limit, c("0", "x1", "X1[0,)"))
      assert_int(x = limit, na.ok = TRUE, lower = 0, null.ok = TRUE)

      # number of tests
      n <- length(private$p_values)

      # designations of tests (if present)
      names <- if(is.data.frame(private$inputs$observations))
        rownames(private$inputs$observations) else
          names(private$inputs$observations)

      if(inputs || pvalue_details || supports){
        pars      <- private$inputs
        idx_alt   <- which(names(pars$computation) == "alternative")
        idx_ex    <- which(names(pars$computation) == "exact")
        idx_dist  <- which(names(pars$computation) == "distribution")
        idx_add   <- which(startsWith(names(pars$computation), "distribution."))
        #idx_cont <- which(names(pars$computation) == "correct")
        idx       <- c(idx_alt, idx_ex, idx_dist, idx_add)#, idx_cont)
      }
      if(supports)
        supp <- self$get_pvalue_supports(unique = FALSE)

      cli_h1(private$test_name)
      cli_text("\n")
      cli_ul(id = "results")
      cli_li("Data: {.field {private$data_name}}")
      if(n > 1) cli_li(paste0("number of tests: {.field {n}}"))

      # determine number of tests to print
      if(is.null(test_idx)) {
        limit <- ifelse(is.null(limit) || is.na(limit), n, limit)
        nums <- seq_len(min(limit, n))
      } else nums <- unique(pmin(na.omit(test_idx), n))

      # number of results to be printed
      len <- length(nums)

      for(i in nums) {
        if(n > 1) {
          cli_end("results")
          cli_text("\n")
          heading <- paste("Results of test", i)
          #cat("Test", i)
          if(!is.null(names) && names[i] != i) {
            # remove line breaks from name
            name <- gsub("[\r\n]", "_", names[i])
            # how many characters can be printed to console (minus reserved)
            chars_tst_line <- getOption("width") - nchar(heading) - 6L
            # how many characters can be used by tag (= row name)
            chars_name_max <- chars_tst_line - 3L
            # length of current row name
            chars_name <- nchar(name)
            # print row name in brackets behind test number
            heading <- paste0(
              heading,
              " (",
              ifelse(
                chars_name > chars_name_max,
                paste(substr(name, 1, chars_name_max - 3L), "..."),
                name
              ),
              ")"
            )
          }
          cli_h2(heading)
          cli_ul(id = "results")
        }

        #cli_ul(id = "test")

        if(inputs) {
          if(is.data.frame(pars$observations)) {
            # fixed number of values (depending on test, e.g. 4 for Fisher's)
            cli_li("{qty(names(pars$observations))} Observation{?s}:")
            # start sub-list
            cli_ul(id = "obs")
            # print vector of strings with names of components and their values
            cli_li(paste0(
              "{.field ",
              names(pars$observations),
              "}: {.val {as.numeric(format(",
              pars$observations[i, ],
              "))}}"
            ))
          } else {
            # observations consist of one or more samples
            cli_li("Observations:")
            # start sub-list
            cli_ul(id = "obs")
            # number of samples per observation
            number_samples <- ifelse(
              is.list(pars$observations[[1]]),
              length(pars$observations),
              1
            )
            # print information about sample(s)
            for(s in seq_len(number_samples)) {
              # current sample
              sample_obs <- if(number_samples > 1)
                pars$observations[[s]][[i]] else
                  pars$observations[[i]]
              # number of observations
              len_obs <- length(sample_obs)
              # number of sample (if more than one) and (some) values
              cli_li(paste(
                "{.field Sample} {ifelse(number_samples > 1,",
                "paste0(col_green(s), ':'), 'of')} {.val {len_obs}}",
                "{mode(sample_obs)} {qty(len_obs)} value{?s}",
                "({ansi_collapse(as.numeric(format(sample_obs)), trunc = 5)})"
              ))
            }
          }
          # end sub-list
          cli_end("obs")

          # print statistics (if any)
          if(
            !is.null(private$statistics) &&
            !all(is.na(private$statistics[i, ]))
          ) {
            cli_li("Statistics:")
            # start sub-list
            cli_ul(id = "statistics")
            # print statistics
            for(j in seq_len(ncol(private$statistics)))
              if(!is.na(private$statistics[i, j]))
                cli_li(paste(
                  "{.field {names(private$statistics)[j]}}:",
                  "{.val {as.numeric(format(private$statistics[i, j]))}}"
                ))
            # end sub-list
            cli_end("statistics")
          }

          # print parameters
          if(!is.null(pars$parameters)) {
            idx_par_valid <- which(!is.na(pars$parameters[i, ]))
            cli_li("Parameters:")
            cli_ul(id = "parameters")
            for(j in idx_par_valid) {
              par_val <- pars$parameters[i, j]
              cli_li(paste(
                "{.field {names(pars$parameters)[j]}}:",
                if(is.numeric(par_val))
                  "{.val {as.numeric(format(par_val))}}" else
                    if(is.logical(par_val))
                      paste0(
                        "{col_blue('", ifelse(par_val, "yes", "no"), "')}"
                      ) else
                        "{.val {par_val}}"
              ))
            }
            cli_end("parameters")
          }

          # print null and alternative hypotheses
          if(pars$computation$alternative[i] == "greater") {
            alt <- "greater than"
            null <- "less than or equal to"
          } else if(pars$computation$alternative[i] == "less") {
            alt <- "less than"
            null <- "greater than or equal to"
          } else {
            alt <- "not equal to"
            null <- "equal to"
          }

          cli_li("Hypotheses:")
          cli_ul(id = "hypotheses")
          cli_li(paste(
            "null: true {names(pars$nullvalues)[1]} is {null}",
            "{.val {pars$nullvalues[[1]][i]}}"
          ))
          cli_li(paste(
            "alternative: true {names(pars$nullvalues)[1]} is {alt}",
            "{.val {pars$nullvalues[[1]][i]}}"
          ))
          cli_end("hypotheses")
        }

        # print p-values
        if(pvalue_details) {
          # print p-values with computation details
          meth <- switch(
            EXPR    = pars$computation$alternative[i],
            minlike = "two-sided, minimum likelihood",
            blaker  = "two-sided, combined tails",
            absdist = "two-sided, absolute distance from mean",
            central = "two-sided, minimum tail doubling",
            less    = "one-sided, lower tail",
            greater = "one-sided, upper tail",
            "two-sided"
          )
          cli_li("p-Value:")
          cli_ul(id = "pv")
          cli_li("computation: {col_blue(meth)}")

          cli_ul(id = "computation")
          cli_li(paste(
           "exact: {ifelse(pars$computation$exact[i],",
           "col_green('yes'), col_red('no'))}"
          ))
          cli_li("distribution: {col_blue(pars$computation$distribution[i])}")
          if(length(idx_add) > 0) {
            cli_ul(id = "dist_pars")
            for(j in idx_add) {
              par_name <- sub("distribution.", "", names(pars$computation)[j])
              par_val <- pars$computation[i, j]
              if(!is.na(par_val))
                cli_li(paste(
                  "{.field {par_name}}:",
                  if(is.numeric(par_val))
                    "{.val {as.numeric(format(par_val))}}" else
                      "{.val {par_val}}"
                ))
            }
            cli_end("dist_pars")
          }

          len_comp <- ncol(pars$computation)
          if(len_comp > length(idx)) {
            idx_comp <- setdiff(seq_len(len_comp), idx)
            if(!all(is.na(pars$computation[i, idx_comp]))) {
              for(j in idx_comp) {
                par_val <- pars$computation[i, j]
                if(!is.na(par_val))
                  cli_li(paste(
                    "{names(pars$computation)[j]}:",
                    if(is.numeric(par_val))
                      "{.val {as.numeric(format(par_val))}}" else
                        if(is.logical(par_val))
                          paste0(
                            "{col_blue('", ifelse(par_val, "yes", "no"), "')}"
                          ) else
                            "{.val {par_val}}"
                  ))
              }
            }
          }
          # if(!pars$computation$exact[i] && !is.na(pars$computation$correct[i]))
          #   cli_li(paste(
          #     "continuity correction: {ifelse(pars$computation$correct[i],",
          #     "col_blue('yes'), col_blue('no'))}"
          #   ))
          cli_end("computation")
          cli_li("value: {.val {as.numeric(format(private$p_values[i]))}}")
          #cli_end("pv")
        } else {
          # print p-values only
          if(supports) {
            cli_li("p-Value:")
            cli_ul(id = "pv")
            cli_li("value: {.val {as.numeric(format(private$p_values[i]))}}")
          } else
            cli_li("p-Value: {.val {as.numeric(format(private$p_values[i]))}}")
        }

        if(supports) {
          cli_li(paste(
            "support: {.val {length(supp[[i]])}} observable p-value{?s}",
            "({ansi_collapse(as.numeric(format(supp[[i]])), trunc = 5)})"
          ))
        }

        if(pvalue_details || supports) cli_end("pv")

        #cli_end("test")
        cli_end()
      }

      if(is.null(test_idx) && limit > 0 && limit < n) {
        cli_text("\n")
        cli_text("\n")

        cli_text(
          paste(
            "[ print limit reached --", n - limit, "results omitted --",
            "use print parameter 'limit' for more results ]\n"
          )
        )
      }
      if(!is.null(test_idx)) {
        cli_text("\n")
        cli_text("\n")

        cli_text(paste("[", length(nums), "out of", n, "results printed ]\n"))
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
