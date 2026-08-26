#' @name wilcox_test_pv
#'
#' @title
#' Wilcoxon's signed-rank test
#'
#' @description
#' `wilcox_test_pv()` performs an exact or approximate Wilcoxon signed-rank test
#' about the location of a population on a single sample or the differences
#' between two paired groups when the data is not necessarily normally
#' distributed. In contrast to [`stats::wilcox.test()`], it is vectorised and
#' only calculates *p*-values. Furthermore, it is capable of returning the
#' discrete *p*-value supports, i.e. all observable *p*-values under a null
#' hypothesis. Multiple tests can be evaluated simultaneously.
#'
#' @param x               numerical vector forming the sample to be tested or a
#'                        list of numerical vectors for multiple tests.
#' @param y               numerical vector forming the second sample to be
#'                        tested or a list of numerical vectors for multiple
#'                        tests; if `y = NULL` (the default), the one-sample
#'                        version is performed; for two-sample tests, all sample
#'                        pairs must have the same length.
#' @param mu              numerical vector or single number of hypothesised
#'                        location(s) for one-sample tests or location shift(s)
#'                        for two-sample tests.
#' @param zero_method     character vector or single string specifying how zero
#'                        differences are handled; must either be `"pratt"` (the
#'                        default) or `"wilcoxon"`. With `"pratt"`, zero
#'                        differences are included when computing ranks, but
#'                        excluded from the test statistic; with `"wilcoxon"`,
#'                        zero differences are discarded entirely before ranking
#'                        (the original Wilcoxon approach).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact_w TRUE
#' @templateVar correct_w TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' The parameters `x`, `mu`, `alternative` and `zero_method` are vectorised. If
#' `x` is a list, they are replicated automatically to have the same lengths. In
#'  case `x` is not a list, it is added to one, which is then replicated to the
#' appropriate length. This allows multiple hypotheses to be tested
#' simultaneously.
#'
#' By setting `exact = NULL`, exact computation is performed only if the sample
#' size in a test setting is smaller than or equal to 200. Otherwise,
#' \eqn{p}-values are computed by normal approximation.
#'
#' The used test statistics `W` is also known as \eqn{T+} and is defined as the
#' sum of ranks of all strictly positive values of the sample `x`.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`], [`mann_whitney_test_pv()`]
#'
#' @references
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. pp. 40-55. \doi{10.1002/9781119196037}
#'
#' @examples
#' # Constructing
#' set.seed(1)
#' r1 <- rnorm(200)
#' r2 <- rnorm(200, 1)
#' r3 <- rnorm(200, 2)
#'
#' ## One-sample tests
#' #  Exact two-sided p-values and their supports
#' results_ex_1s <- wilcox_test_pv(r1)
#' print(results_ex_1s)
#' results_ex_1s$get_pvalues()
#' results_ex_1s$get_pvalue_supports()
#'
#' #  Multiple normal-approximated one-sided tests ("greater")
#' results_ap_1s <- wilcox_test_pv(list(r1, r2), alternative = "greater", exact = FALSE)
#' print(results_ap_1s)
#' results_ap_1s$get_pvalues()
#' results_ap_1s$get_pvalue_supports()
#'
#' ## Two-sample-tests
#' #  Normal-approximated one-sided p-values ("less") and their supports
#' results_ap_2s <- wilcox_test_pv(r1, r2, alternative = "less", exact = FALSE)
#' print(results_ap_2s)
#' results_ap_2s$get_pvalues()
#' results_ap_2s$get_pvalue_supports()
#'
#' #  Multiple exact two-sided tests ("greater")
#' results_ex_2s <- wilcox_test_pv(list(r1, r3), r2)
#' print(results_ex_2s)
#' results_ex_2s$get_pvalues()
#' results_ex_2s$get_pvalue_supports()
#'
#' @importFrom checkmate qassert qassertr
#' @importFrom cli cli_abort cli_warn
#' @export
wilcox_test_pv <- function(
  x,
  y = NULL,
  mu = 0,
  alternative = "two.sided",
  exact = NULL,
  correct = TRUE,
  digits_rank = Inf,
  zero_method = "pratt",
  simple_output = FALSE
) {
  # plausibility checks of input parameters
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  one_sample <- is.null(y)
  if(!one_sample) {
    qassert(y, c("N+", "L+"))
    if(!is.list(y)) y <- list(y) else qassertr(y, "N+")
  }
  len_y <- length(y)

  qassert(mu, "N+()")
  len_m <- length(mu)

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      tolower(alternative[i]),
      c("two.sided", "less", "greater")
    )
  }

  qassert(exact,   c("B1", "0"))

  qassert(correct, c("B1", "X1[0, 3]"))
  if(!is.logical(correct) && (is.null(exact) || (!is.null(exact) && !exact))) {
    edgeworth <- round(correct)
    correct   <- TRUE
  } else {
    if(!is.null(exact) && exact && is.numeric(correct)) correct <- TRUE
    edgeworth <- 0
  }

  qassert(digits_rank, "N1")

  len_z <- length(zero_method)
  for(i in seq_len(len_z)){
    zero_method[i] <- match.arg(
      tolower(zero_method[i]),
      c("pratt", "wilcoxon")
    )
  }

  qassert(simple_output, "B1")

  # replicate inputs to same length
  len_g                         <- max(len_x, len_y, len_m, len_a, len_z)
  if(len_x < len_g) x           <- rep_len(x, len_g)
  if(len_m < len_g) mu          <- rep_len(mu, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)
  if(len_z < len_g) zero_method <- rep_len(zero_method, len_g)

  if(!one_sample) {
    if(len_y < len_g) y <- rep_len(y, len_g)
    res_obs <- list(x, y)

    for(i in seq_len(len_g)) {
      # check if lengths of all sample pairs are equal; stop if they are not
      if(length(x[[i]]) != length(y[[i]]))
        cli_abort('All paired samples must have the same length')
      # compute differences
      x[[i]] <- x[[i]] - y[[i]]
    }
  }

  # compute ranks and lengths
  n     <- integer(len_g)
  m     <- integer(len_g)
  W     <- numeric(len_g)
  means <- numeric(len_g)
  sds   <- numeric(len_g)
  ties  <- logical(len_g)
  zeros <- numeric(len_g)
  stats <- vector("list", len_g)
  for(i in seq_len(len_g)) {
    y <- x[[i]] - mu[i]
    m[i] <- length(y)

    is_zero  <- (y == 0)
    zeros[i] <- sum(is_zero)
    if(zero_method[i] == "wilcoxon") {
      y         <- y[!is_zero]
      is_zero   <- rep(FALSE, length(y))
      n[i]      <- m[i] - zeros[i]
    } else n[i] <- m[i]

    ranks <- if(is.finite(digits_rank))
      rank(abs(signif(y, digits_rank))) else
        rank(abs(y))

    ranks   <- ranks[!is_zero]
    pos_y   <- y[!is_zero] > 0
    ties[i] <- m[i] - zeros[i] != length(unique(ranks))

    # test statistics
    W[i] <- sum(ranks[pos_y])
    if(ties[i] || (zeros[i] & zero_method[i] == "pratt") ||
       is.null(exact) && n > 200 || !is.null(exact) && !exact)
      stats[[i]] <- as.integer(round(2 * ranks))

    # parameters for normal approximation
    means[i] <- n[i] * (n[i] + 1) / 4
    if(zeros[i] && zero_method[i] == "pratt")
      means[i] <- means[i] - zeros[i] * (zeros[i] + 1) / 4
    if(is.null(exact) && n[i] > 200 || !is.null(exact) && !exact) {
      sds[i] <- means[i] * (2 * n[i] + 1) / 6
      # correct for zeros depending on `zero_method`
      if(zeros[i] && zero_method[i] == "pratt") {
        sds[i] <- sds[i] - zeros[i] * (zeros[i] + 1) * (2 * zeros[i] + 1) / 24
      }
      # correct for possible ties
      if(ties[i]) {
        t <- table(ranks)
        sds[i] <- sds[i] - sum(t^3 - t) / 48
      }
    }
  }
  sds <- sqrt(sds)

  ex <- if(is.null(exact)) n <= 200 else rep(exact, len_g)
  ew <- edgeworth & !ex
  ew[ex] <- NA

  # compute Edgeworth coefficients for normal approximations, if desired
  idx_ew <- which(ew)
  if(length(idx_ew)) {
    ew_coefs <- matrix(NA_real_, len_g, edgeworth)
    for(i in idx_ew) {
      if(edgeworth >= 1) {
        w <- stats[[i]]/2
        s <- sqrt(sum(w^2) / 4)
        ew_coefs[i, 1] <- -sum(w^4) / s^4 / 192
      }
      if(edgeworth >= 2)
        ew_coefs[i, 2] <- sum(w^6) / s^6 / 2880
      if(edgeworth == 3)
        ew_coefs[i, 3] <- ew_coefs[i, 1]^2 / 2
    }
  }

  # determine unique parameter sets
  params    <- data.frame(
    alternative, n, means, sds, ew, ties, zeros, zero_method
  )
  params_ex <- unique(
    subset(params, ex & !ties & (!zeros | zero_method == "wilcoxon"), 1:2)
  )
  params_tz <- subset(params, ex & (ties | zeros & zero_method == "pratt"))
  params_ap <- unique(subset(params, !ex, -2))
  idx_ex    <- as.numeric(rownames(params_ex))
  idx_tz    <- as.numeric(rownames(params_tz))
  idx_ap    <- as.numeric(rownames(params_ap))
  rows      <- c(idx_ex, idx_tz, idx_ap)
  params_u  <- params[rows, ]
  if(any(ew, na.rm = TRUE)) ew_coefs_u <- ew_coefs[rows, , drop = FALSE]

  len_ex <- length(idx_ex)
  len_tz <- length(idx_tz)
  len_ap <- length(idx_ap)
  idx_ex <- seq_len(len_ex)
  idx_ap <- len_ex + len_tz + seq_len(len_ap)
  len_u  <- len_ex + len_tz + len_ap

  alts_u      <- params_u$alternative
  n_u         <- params_u$n
  mean_u      <- params_u$means
  sd_u        <- params_u$sds
  ew_u        <- params_u$ew
  ties_u      <- params_u$ties
  zeros_u     <- params_u$zeros
  zero_meth_u <- params_u$zero_method

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # pre-compute exact distributions (if any)
  sizes_ex <- unique(n_u[idx_ex])
  d <- generate_sign_rank_probs(sizes_ex)

  # begin exact computations (if any)
  for(i in idx_ex) {
    idx_supp <- which(
      alts_u[i] == alternative & n_u[i] == n &
        ex & !ties & (!zeros | zero_method == "wilcoxon")
    )

    idx_d <- which(sizes_ex == n_u[i])
    d_i <- numerical_adjust(d[[idx_d]])

    if(simple_output) {
      # compute p-values directly
      res[idx_supp] <- switch(
        EXPR = alts_u[i],
        less = p_from_d(W[idx_supp], d_i),
        greater = p_from_d(W[idx_supp] - 1, d_i, FALSE),
        two.sided = {
          idx_l <- which(W[idx_supp] < mean_u[i])
          idx_u <- which(W[idx_supp] >= mean_u[i])
          pv <- numeric(length(idx_supp))
          if(length(idx_l))
            pv[idx_l] <- p_from_d(W[idx_supp][idx_l], d_i)
          if(length(idx_u))
            pv[idx_u] <- p_from_d(W[idx_supp][idx_u] - 1, d_i, FALSE)
          pmin(1, 2 * pv)
        }
      )
    } else {
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[i],
        probs = d_i,
        expectations = abs(seq_along(d) - 1 - mean_u[i])
      )

      # store results and support
      res[idx_supp] <- pv_supp[W[idx_supp] + 1]
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  # begin exact computations for settings with ties or zeros (if any)
  idx_out <- len_ex + 1
  for(i in idx_tz) {
    # compute exact distribution with ties
    dist <- numerical_adjust(sign_rank_probs_ties_int(stats[[i]]))
    R <- 2 * W[i]

    if(simple_output) {
      # compute p-values directly
      res[i] <- switch(
        EXPR = alternative[i],
        less = p_from_d(R, dist),
        greater = p_from_d(R - 1, dist, FALSE),
        two.sided = pmin(1, 2 * p_from_d(
          R - W[i] > means[i], dist, W[i] <= means[i]
        ))
      )
    } else {
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alternative[i],
        probs = dist
      )

      # store results and support
      res[i]              <- pv_supp[R + 1]
      supports[[idx_out]] <- unique(sort(pv_supp))
      indices[[idx_out]]  <- i
      idx_out             <- idx_out + 1
    }
  }

  # begin approximation computations (if any)
  for(i in idx_ap) {
    idx_supp <- which(
      alts_u[i] == alternative &
      !ex &
      ties_u[i] == ties &
      zeros_u[i] == zeros &
      zero_meth_u[i] == zero_method &
      mean_u[i] == means &
      sd_u[i] == sds
    )

    ew_ok <- !is.na(ew_u[i]) && ew_u[i]
    e <- if(ew_ok) ew_coefs_u[i, ]

    if(simple_output) {
      res[idx_supp] <- switch(
        EXPR = alts_u[i],
        less = pnorm_MW_edgeworth(
          W[idx_supp], mean_u[i], sd_u[i], TRUE, correct, e
        ),
        greater = pnorm_MW_edgeworth(
          W[idx_supp], mean_u[i], sd_u[i], FALSE, correct, e
        ),
        two.sided = pmin(1, 2 * if(ew_ok)
          pmin(
            pnorm_MW_edgeworth(
              W[idx_supp], mean_u[i], sd_u[i], TRUE, correct, e
            ),
            pnorm_MW_edgeworth(
              W[idx_supp], mean_u[i], sd_u[i], FALSE, correct, e
            )
          ) else pmin(1,
            pnorm(-abs(W[idx_supp] - mean_u[i]), -correct * 0.5, sd_u[i])
          )
        )
      )
    } else {
      # compute p-value support
      z <- if(!any(ties[idx_supp])) 0L:(n_u[i] * (n_u[i] + 1L)) else
        0L:(n_u[i] * (n_u[i] + 1L)) / 2
      pv_supp <- switch(
        EXPR = alts_u[i],
        less = pnorm_MW_edgeworth(z, mean_u[i], sd_u[i], TRUE, correct, e),
        greater = pnorm_MW_edgeworth(z, mean_u[i], sd_u[i], FALSE, correct, e),
        two.sided = pmin(1, 2 * if(ew_ok)
          pmin(
            pnorm_MW_edgeworth(z, mean_u[i], sd_u[i], TRUE, correct, e),
            pnorm_MW_edgeworth(z, mean_u[i], sd_u[i], FALSE, correct, e)
          ) else pnorm(-abs(z - mean_u[i]), -correct * 0.5, sd_u[i])
        )
      )

      # store results and support
      idx_stat <- 1 + if(!any(ties[idx_supp]))
        W[idx_supp] else round(2 * W[idx_supp])
      res[idx_supp] <- pv_supp[idx_stat]
      if(!simple_output) {
        supports[[i]] <- unique(sort(pv_supp))
        indices[[i]]  <- idx_supp
      }
    }
  }

  # remove empty supports and indices
  if(!simple_output) {
    ord <- order(sapply(indices, "[", 1))
    supports <- supports[ord]
    indices  <- indices[ord]
  }

  # create output object
  out <- if(!simple_output) {
    arg_names  <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon signed-rank test",
      inputs = list(
        observations = if(one_sample) x else res_obs,
        parameters = NULL,
        nullvalues = if(one_sample) data.frame(location = mu) else
          data.frame(`location shift` = mu, check.names = FALSE),
        computation = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            alternative = alternative,
            exact = ex,
            distribution = ifelse(
              ex,
              paste0(
                "Wilcoxon's signed-rank",
                ifelse(ties | zeros & zero_method == "pratt", " (", ""),
                ifelse(ties, "tie-", ""),
                ifelse(ties & zeros & zero_method == "pratt", " and ", ""),
                ifelse(zeros & zero_method == "pratt", "zero-", ""),
                ifelse(ties | zeros & zero_method == "pratt", "adjusted)", "")
              ),
              "normal"
            ),
            #distribution.mean = ifelse(!ex, means, NA_real_),
            #distribution.sd = ifelse(!ex, sds, NA_real_),
            `continuity correction` = ifelse(ex, NA, correct),
            `Edgeworth expansion` = ew,
            `Edgeworth series terms` = ifelse(ew, ew * edgeworth, NA),
            `sample size` = m,
            zeros = ifelse(zeros, zeros, FALSE),
            `zero handling` = ifelse(
              zero_method == "pratt", "Pratt", "Wilcoxon"
            ),
            `(non-zero) ties` = ties,
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(W),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = arg_names[if(one_sample) "x" else c("x", "y")]
    )
  } else res

  # return results
  return(out)
}
