#' @name perm_test_pv
#'
#' @title
#' Permutation Test
#'
#' @description
#' `perm_test_pv()` performs a permutation test for the comparison of two
#' independent samples using a user-supplied test statistic. It is vectorised
#' over multiple test problems: by passing lists of sample vectors for `x` and
#' `y`, several tests can be evaluated simultaneously. The function can compute
#' exact *p*-values (by enumerating all permutations) or approximate *p*-values
#' via Monte Carlo sampling. It is capable of returning the discrete *p*-value
#' supports, i.e. all observable *p*-values under the null hypothesis.
#'
#' @param x                 numerical vector or list of numerical vectors of
#'                          observations in the first group.
#' @param y                 numerical vector or list of numerical vectors of
#'                          observations in the second group.
#' @param statistic         single character giving the type of test statistic
#'                          to be used; must be one of `"diff_mean"`,
#'                          `"diff_median"`, `"diff_t"`, `"diff_hl"`,
#'                          `"ratio_var"` or `"ratio_sd"` (see Details).
#' @param mu                numerical vector giving the hypothesised values
#'                          under the null; if `mu == NULL` (the default) it
#'                          is set to 0 for test statistics that are based on
#'                          differences (`diff_`*) or 1 for ratio-based tests
#'                          statistics (`ratio_`*).
#' @param exact             logical value that indicates whether \eqn{p}-values
#'                          are to be calculated by exact computation (`TRUE`)
#'                          or by simulation (`FALSE`); if `NULL` (the default)
#'                          exact computation is performed if the number of
#'                          possible sample combinations (see Details) is below
#'                          or equal to `max_exact_combs`; otherwise the
#'                          respective \eqn{p}-values are obtained by
#'                          simulation.
#' @param max_exact_combs   maximum number of allowed combinations for exact
#'                          computation of the permutation distribution (see
#'                          details).
#' @param MC_sims           positive integer or `NULL`. The number of Monte
#'                          Carlo permutations to draw when `exact = FALSE`.
#'                          Ignored when `exact = TRUE`. Default is `10000L`.
#' @param seed              single integer or `NULL`. Random seed passed to
#'                          [`set.seed()`] before Monte Carlo sampling to ensure
#'                          reproducibility. Use `NULL` to skip seed setting.
#'                          Ignored when `exact = TRUE`.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar simple.output TRUE
#'
#' @details
#' Under the null hypothesis of exchangeability (i.e. both samples come from
#' the same distribution), all \eqn{\binom{n_x + n_y}{n_x}} partitions of the
#' pooled sample into groups of sizes \eqn{n_x} and \eqn{n_y} are equally
#' likely. `perm_test_pv()` exploits this by:
#'
#' \enumerate{
#'   \item Computing the observed test statistic
#'     \eqn{T_\text{obs} = \texttt{stat_fun}(x, y)}.
#'   \item Enumerating (if `exact = TRUE`) or randomly sampling (if
#'     `exact = FALSE`) permutations of the pooled sample.
#'   \item Evaluating `stat_fun()` on each permuted split.
#'   \item Deriving the *p*-value as the fraction of permutation statistics
#'     that are at least as extreme as \eqn{T_\text{obs}}, where *extreme*
#'     is defined by `alternative`:
#'     \describe{
#'       \item{`"less"`}{fraction with permuted statistic
#'         \eqn{\le T_\text{obs}}}
#'       \item{`"greater"`}{fraction with permuted statistic
#'         \eqn{\ge T_\text{obs}}}
#'       \item{`"two.sided"`}{fraction with
#'         \eqn{|\text{permuted statistic}| \ge |T_\text{obs}|}}
#'     }
#' }
#'
#' Exact computation is only feasible for small pooled samples; the total
#' number of permutations grows as \eqn{\binom{n_x + n_y}{n_x}}. If
#' `exact = NULL`, exact computation is applied when the number of computations
#' is below or equal to `max_exact_combs` (default \eqn{20{,}000}). Otherwise,
#' the distribution of the test statistics is approximated by Monte Carlo
#' simulation. When `exact = TRUE` and the number of permutations exceeds
#' `max_exact_combs`, an error is raised. Set `exact = FALSE` to use Monte Carlo
#' sampling in such cases.
#'
#' The parameters `x`, `y`, `alternative`, and `MC_sims` are vectorised via
#' lists: if `x` and `y` are lists of the same length, each pair
#' `(x[[i]], y[[i]])` defines one test. Scalars and single vectors are
#' recycled to the required length.
#'
#' @template return
#'
#' @seealso
#' [`mann_whitney_test_pv()`], [`stats::wilcox.test()`]
#'
#' @references
#' Pitman, E. J. G. (1937). Significance tests which may be applied to
#'   samples from any populations. *Supplement to the Journal of the Royal
#'   Statistical Society*, **4**(1), pp. 119–130. \doi{10.2307/2984124}
#'
#' Good, P. (2000). *Permutation Tests: A Practical Guide to Resampling
#'   Methods for Testing Hypotheses*. Second Edition. New York: Springer.
#'   \doi{10.1007/978-1-4757-3235-1}
#'
#' @examples
#' set.seed(42)
#' x1 <- rnorm(8)
#' y1 <- rnorm(8, mean = 1)
#'
#' # Exact two-sided p-value (default: difference of means)
#' results_ex <- perm_test_pv(x1, y1)
#' print(results_ex)
#' results_ex$get_pvalues()
#' results_ex$get_pvalue_supports()
#'
#' # Monte Carlo approximation with a custom statistic (median difference)
#' results_mc <- perm_test_pv(
#'   x1, y1,
#'   stat_fun = function(x, y) median(x) - median(y),
#'   exact = FALSE, MC_sims = 9999L, seed = 1L
#' )
#' results_mc$get_pvalues()
#'
#' # Multiple tests simultaneously
#' xs <- list(rnorm(6), rnorm(5, 1))
#' ys <- list(rnorm(6, 1), rnorm(5, 2))
#' results_v <- perm_test_pv(xs, ys, exact = FALSE, MC_sims = 4999L, seed = 7L)
#' results_v$get_pvalues()
#'
#' @importFrom checkmate qassert qassertr
#' @importFrom cli cli_warn
#' @export
perm_test_pv <- function(
  x,
  y,
  statistic = c(
    "diff_mean", "diff_median", "diff_t", "diff_hl", "ratio_var", "ratio_sd"
  ),
  mu = NULL,
  alternative = "two.sided",
  exact = NULL,
  max_exact_combs = 10L^7,
  MC_sims = 10000L,
  seed = NULL,
  simple_output = FALSE
) {
  # input checks
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  qassert(y, c("N+", "L+"))
  if(!is.list(y)) y <- list(y) else qassertr(y, "N+")
  len_y <- length(y)

  statistic <- match.arg(
    tolower(statistic),
    c("diff_mean", "diff_median", "diff_t", "diff_hl", "ratio_var", "ratio_sd")
  )

  if(startsWith(statistic, "diff")) {
    qassert(mu, c("n+()", "0"))
    if(!is.null(mu)) mu[is.na(mu)] <- 0 else mu <- 0
  } else {
    qassert(mu, c("n+(0,)", "0"))
    if(!is.null(mu)) mu[is.na(mu)] <- 1 else mu <- 1
  }
  len_m <- length(mu)

  len_a <- length(alternative)
  for(i in seq_len(len_a)) {
    alternative[i] <- match.arg(
      tolower(alternative[i]),
      c("two.sided", "less", "greater")
    )
  }

  qassert(exact, c("B1", "0"))

  qassert(max_exact_combs, "X1[1,)")
  max_exact_combs <- as.integer(max_exact_combs)

  qassert(MC_sims, "X1[1,)")
  MC_sims <- as.integer(MC_sims)

  assert_numeric(seed, any.missing = FALSE, null.ok = TRUE)

  qassert(simple_output, "B1")

  # replicate to common length
  len_g <- max(len_x, len_y, len_m, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_y < len_g) y <- rep_len(y, len_g)
  if(len_m < len_g) mu <- rep_len(mu, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # set test statistics function according to 'statistic' parameter
  stat_fun <- switch(
    statistic,
    diff_mean = function(x, y) mean(x) - mean(y),
    diff_median = function(x, y) median(x) - median(y),
    diff_t = function(x, y) {
      mx <- length(x)
      my <- length(y)
      sp <- sqrt(((mx - 1) * var(x) + (my - 1) * var(y)) / (mx + my - 2))
      (mean(x) - mean(y)) / sp / sqrt(1/mx + 1/my)
    },
    diff_hl = function(x, y) median(outer(x, y, "-")),
    ratio_var = function(x, y) var(x) / var(y),
    ratio_sd = function(x, y) sd(x) / sd(y)
  )

  # per-test sizes and observed statistics
  nx <- integer(len_g)
  ny <- integer(len_g)
  T_obs <- numeric(len_g)
  n_total <- integer(len_g)
  use_exact <- logical(len_g)

  for(i in seq_len(len_g)) {
    nx[i] <- length(x[[i]])
    ny[i] <- length(y[[i]])
    n_total[i] <- nx[i] + ny[i]

    # check whether exact computation is feasible
    n_perm_i <- choose(n_total[i], nx[i])
    use_exact[i] <- if(is.null(exact)) n_perm_i <= max_exact_combs else exact

    if(is.null(exact) && n_perm_i > max_exact_combs)
      cli_warn(
        c(
          paste0(
            "Test ", i, ": exact computation requires ",
            formatC(n_perm_i, format = "fg", big.mark = ","),
            " permutations, exceeding the limit given by 'max_exact_combs' of ",
            formatC(max_exact_combs, format = "fg", big.mark = ","), "."
          ),
          i = paste(
            "Increase limit for exact calculation.",
            "Using Monte Carlo sampling instead."
          )
        )
      )

    # modify 'x' according to null
    x_null <- switch(
      statistic,
      ratio_var = x[[i]] / sqrt(mu[i]),
      ratio_sd  = x[[i]] / mu[i],
      x[[i]] - mu[i]
    )

    # observed statistic
    T_obs[i] <- stat_fun(x_null, y[[i]])
  }

  if(!is.null(seed) && any(!use_exact)) set.seed(seed)

  # prepare output containers
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_g)
    indices  <- vector("list", len_g)
  }

  # compute permutation p-values
  for(i in seq_len(len_g)) {
    if(use_exact[i]) {
      # exact: enumerate all splits
      exact_comp <- perm_test_all_combs(x[[i]], y[[i]], mu[i], statistic)

      # observable statistics (unique and sorted!)
      T_stats <- exact_comp$statistics

      # adjust observable statistics for slight numerical differences
      T_stats <- numerical_adjust(T_stats, FALSE)

      # observed statistics
      T_obs[i] <- T_stats[which.min(abs(T_stats - T_obs[i]))]#exact_comp$observed

      # adjust test statistics if i-th test is two-sided
      if(alternative[i] == "two.sided") {
        T_stats <- abs(T_stats)
        T_obs[i] <- abs(T_obs[i])
      }

      # probabilities of each statistic
      probs <- exact_comp$probabilities

      # order statistics (and probabilities accordingly)
      ord_T   <- order(T_stats)
      T_stats <- T_stats[ord_T]
      probs   <- probs[ord_T]

      # adjust probabilities accordingly (if necessary)
      dupl <- duplicated(T_stats)
      if(any(dupl)) {
        sf <- stepfun(T_stats, c(0, cumsum(probs)))
        T_stats <- T_stats[!dupl]
        probs <- diff(c(0, sf(T_stats)))
      }

      # adjust probabilities for slight numerical differences
      probs <- numerical_adjust(probs)
    } else {
      # pooled observations
      pooled <- c(x_null, y[[i]])

      # Monte Carlo: draw MC_sims random permutations
      T_stats <- c(rep(0, MC_sims - 1), T_obs[i])
      for(j in seq_len(MC_sims - 1)) {
        perm <- sample.int(n_total[i])
        xi <- pooled[perm[seq_len(nx[i])]]
        yi <- pooled[perm[seq_len(ny[i]) + nx[i]]]
        T_stats[j] <- stat_fun(xi, yi)
      }

      # adjust test statistics if i-th test is two-sided
      if(alternative[i] == "two.sided") T_stats <- abs(T_stats)

      # adjust observable statistics for slight numerical differences
      T_stats <- numerical_adjust(T_stats, FALSE)
      T_obs[i] <- T_stats[MC_sims]

      # compute probabilities and adjust for slight numerical differences
      probs <- numerical_adjust(as.numeric(prop.table(table(T_stats))))

      # sorted and unique statistics
      T_stats <- unique(sort(T_stats))
    }

    # unique sorted p-value support from all permutation statistics
    pv_supp <- switch(
      alternative[i],
      less = cumsum(probs),                 #sapply(T_perm_sorted, function(t) mean(T_stats <= t)),
      greater = rev(cumsum(rev(probs))),    #sapply(T_perm_sorted, function(t) mean(T_stats >= t)),
      two.sided = rev(cumsum(rev(probs)))   #sapply(T_perm_sorted, function(t) mean(abs(T_stats) >= abs(t)))
    )

    idx_pv <- which(T_stats == T_obs[i])
    res[i] <- pv_supp[idx_pv]
    if(!simple_output) {
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- i
    }
  }

  # create output object
  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)
    stat_name <- if(is.na(dnames["statistic"]))
      deparse1(substitute(statistic)) else dnames["statistic"]
    switch(
      statistic,
      diff_mean = {
        null_label <- "location shift"
        stat_label <- "difference of means"
      },
      diff_median = {
        null_label <- "location shift"
        stat_label <- "difference of medians"
      },
      diff_t = {
        null_label <- "location shift"
        stat_label <- "Student (t)"
      },
      diff_hl = {
        null_label <- "location shift"
        stat_label <- "Hodges-Lehmann"
      },
      ratio_var = {
        null_label <- "ratio of variances"
        stat_label <- "ratio of variances"
      },
      ratio_sd = {
        null_label <- "ratio of standard deviations"
        stat_label <- "ratio of standard deviations"
      }
    )

    DiscreteTestResults$new(
      test_name = "Permutation test",
      inputs = list(
        observations = list(x, y),
        parameters   = data.frame(
          `size of first sample` = nx,
          `size of second sample` = ny,
          check.names = FALSE
        ),
        nullvalues = data.frame(
          setNames(list(mu), null_label),
          check.names = FALSE
        ),
        computation  = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            alternative = alternative,
            exact = use_exact,
            distribution = ifelse(
              use_exact,
              paste(
                "permutation (exact,",
                formatC(choose(n_total, nx[i]), format = "fg", big.mark = ","),
                "combinations)"
              ),
              paste0("permutation (MC, B = ", MC_sims, ")")),
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(
        `observed statistic` = T_obs,
        `test statistic type` = stat_label,
        check.names = FALSE
      ),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = dnames[c("x", "y")]
    )
  } else res

  # return results
  return(out)
}
