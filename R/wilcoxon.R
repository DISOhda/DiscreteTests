#' @name wilcoxon_test_pv
#'
#' @title
#' Wilcoxon sign rank test
#'
#' @description
#' `wilcoxon_test_pv()` performs an exact or approximate Wilcoxon sign rank test
#' about the location of a population on a single sample. In contrast to
#' [`stats::wilcox.test()`], only the one-sample test is performed, but it is
#' vectorised and only calculates *p*-values. Furthermore, it is capable of
#' returning the discrete *p*-value supports, i.e. all observable *p*-values
#' under a null hypothesis. Multiple tests can be evaluated simultaneously.
#'
#' @param x             numerical vector forming the sample to be tested or list
#'                      of numerical vectors for multiple samples.
#' @param mu            numerical vector of hypothesised location(s).
#' @param digits_rank   single number giving the significant digits used to
#'                      compute ranks for the test statistics.
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#'
#' @details
#' The parameters `x`, `mu` and `alternative` are vectorised. If `x` is a list,
#' they are replicated automatically to have the same lengths. In case `x` is
#' not a list, it is added to one, which is then replicated to the appropriate
#' length. This allows multiple hypotheses to be tested simultaneously.
#'
#' In the presence of ties or observations that are equal to `mu`, computation
#' of exact *p*-values is not possible. Therefore, `exact` is ignored in these
#' cases and *p*-values are calculated by a normal approximation.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`]
#'
#' @references
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. pp. 40-55. \doi{10.1002/9781119196037}
#'
#' @examples
#' # Constructing
#' set.seed(1)
#' r <- rnorm(1000)
#'
#' # Computation of exact two-sided p-values and their supports
#' results_ex  <- wilcox_test_pv(r)
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- wilcox_test_pv(r, alternative = "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom stats psignrank
#' @importFrom checkmate qassert qassertr
#' @export
wilcox_test_pv <- function(
  x,
  mu = 0,
  alternative = "two.sided",
  exact = TRUE,
  correct = TRUE,
  digits_rank = Inf,
  simple_output = FALSE
) {
  # catch input values or data names from call
  dnames <- sapply(match.call(), deparse1)
  if(is.na(dnames["mu"])) dnames <- c(dnames, mu = deparse1(substitute(mu)))

  # plausibility checks of input parameters
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  qassert(mu, "N+()")
  len_m <- length(mu)

  qassert(exact, "B1")
  if(!exact) qassert(correct, "B1")

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
  }

  qassert(digits_rank, "N")

  qassert(simple_output, "B1")

  # replicate inputs to same length
  len_g <- max(len_x, len_m, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_m < len_g) mu <- rep_len(mu, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # compute ranks and lengths
  n <- integer(len_g)
  #ranks <- vector("list", len_g)
  V <- numeric(len_g)
  means <- numeric(len_g)
  vars <- numeric(len_g)
  any_zeros <- logical(len_g)
  any_ties <- logical(len_g)
  for(i in seq_len(len_g)) {
    x[[i]] <- x[[i]] - mu[i]

    is_zero <- (x[[i]] == 0)
    any_zeros[i] <- any(is_zero)
    if(any_zeros[i]) x[[i]] <- x[[i]][!is_zero]

    n[i] <- length(x[[i]])

    ranks <- if(is.finite(digits_rank))
      rank(abs(signif(x[[i]], digits_rank))) else
        rank(abs(x[[i]]))

    V[i] <- sum(ranks[x[[i]] > 0])
    any_ties[i] <- length(ranks) != length(unique(ranks))

    means[i] <- n[i] * (n[i] + 1) / 4
    t <- table(ranks)
    vars[i] <- sqrt(n[i] * (n[i] + 1) * (2 * n[i] + 1) / 24 - sum(t^3 - t) / 48)
  }

  # determine unique parameter sets
  params <- data.frame(n, any_zeros, any_ties, alternative, means, vars)
  if(exact) {
    params_ex <- unique(subset(params, !any_zeros & !any_ties, c(1, 4)))
    params_ap <- unique(subset(params, any_zeros | any_ties, -(1:3)))
  } else {
    params_ex <- params[NULL, c(1, 4)]
    params_ap <- unique(params[-(1:3)])
  }
  idx_ex   <- as.numeric(rownames(params_ex))
  idx_ap   <- as.numeric(rownames(params_ap))
  rows     <- c(idx_ex, idx_ap)
  params_u <- params[rows, ]

  len_ex <- length(idx_ex)
  len_ap <- length(idx_ap)
  idx_ex <- seq_len(len_ex)
  idx_ap <- len_ex + seq_len(len_ap)
  len_u  <- len_ex + len_ap

  n_u     <- params_u$n
  zeros_u <- params_u$any_zeros
  ties_u  <- params_u$any_ties
  alts_u  <- params_u$alternative
  means_u <- params_u$means
  vars_u  <- params_u$vars

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # begin exact computations (if any)
  for(i in idx_ex) {
    idx_supp <- which(n_u[i] == n & alts_u[i] == alternative &
                        !any_zeros & !any_ties)

    if(simple_output) {
      # compute p-values directly
      res[idx_supp] <- switch(
        EXPR = alts_u[i],
        less = psignrank(V[idx_supp], n_u[i]),
        greater = psignrank(V[idx_supp] - 1, n_u[i], lower.tail = FALSE),
        two.sided = {
          idx_l <- which(V[idx_supp] < means_u[i])
          idx_u <- which(V[idx_supp] >= means_u[i])
          pv <- numeric(length(idx_supp))
          if(length(idx_l))
            pv[idx_l] <- psignrank(V[idx_supp][idx_l], n_u[i])
          if(length(idx_u))
            pv[idx_u] <- psignrank(n_u[i] - V[idx_supp][idx_u], n_u[i])
          pmin(1, 2 * pv)
        }
      )
    } else {
      # generate all probabilities under current n
      d <- generate_signrank_probs(n_u[i])
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[i],
        probs = d,
        expectations = abs(seq_along(d) - 1 - means_u[i])
      )

      # store results and support
      res[idx_supp] <- pv_supp[V[idx_supp] + 1]
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  # begin approximation computations (if any)
  for(i in idx_ap) {
    if(exact) {
      idx_supp <- which(alts_u[i] == alternative & (any_zeros | any_ties) &
                         means_u[i] == means & vars_u[i] == vars)
    } else {
      idx_supp <- which(alts_u[i] == alternative &
                          means_u[i] == means & vars_u[i] == vars)
    }

    if(simple_output) {
      res[idx_supp] <- support_normal(
        alternative = alts_u[i],
        x = V[idx_supp],
        mean = means_u[i],
        sd = vars_u[i],
        correct = correct
      )
    } else {
      # compute p-value support
      pv_supp <- support_normal(
        alternative = alts_u[i],
        x = 0L:((n_u[i] * (n_u[i] + 1L)) %/% 2L),
        mean = means_u[i],
        sd = vars_u[i],
        correct = correct
      )

      # store results and support
      res[idx_supp] <- pv_supp[V[idx_supp] + 1]
      if(!simple_output) {
        supports[[i]] <- unique(sort(pv_supp))
        indices[[i]]  <- idx_supp
      }
    }
  }

  out <- if(!simple_output) {
    exact_v <- exact & !any_zeros & !any_ties

    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Wilcoxon signed rank test",
        paste0(
          "Normal-approximated Wilcoxon signed rank test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(
        observations = x,
        nullvalues = data.frame(location = mu),
        parameters = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            n = ifelse(exact_v, n, NA),
            mean = ifelse(!exact_v, means, NA),
            sd = ifelse(!exact_v, sqrt(vars), NA),
            alternative = alternative,
            exact = exact_v,
            distribution = ifelse(exact_v, "Wilcoxon's sign rank", "normal")
          )
        )
      ),
      statistics = data.frame(V),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = paste0(dnames["x"], " and ", dnames["mu"])
    )
  } else res

  return(out)
}
