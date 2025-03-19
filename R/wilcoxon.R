#' @name wilcox_single_test_pv
#'
#' @title
#' Wilcoxon one-sample signed-rank test
#'
#' @description
#' `wilcox_single_test_pv()` performs an exact or approximate Wilcoxon sign rank
#' test about the location of a population on a **single** sample. In contrast
#' to [`stats::wilcox.test()`], only the one-sample test is performed, but it is
#' vectorised and only calculates *p*-values. Furthermore, it is capable of
#' returning the discrete *p*-value supports, i.e. all observable *p*-values
#' under a null hypothesis. Multiple tests can be evaluated simultaneously.
#'
#' @param x             numerical vector forming the sample to be tested or list
#'                      of numerical vectors for multiple samples.
#' @param mu            numerical vector of hypothesised location(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' The parameters `x`, `mu` and `alternative` are vectorised. If `x` is a list,
#' they are replicated automatically to have the same lengths. In case `x` is
#' not a list, it is added to one, which is then replicated to the appropriate
#' length. This allows multiple hypotheses to be tested simultaneously.
#'
#' In the presence of ties or observations that are equal to `mu`, computation
#' of exact *p*-values is not possible. This also applies, if the sample size is
#' greater than 1,038, because [`stats::dsignrank`] then produces `Inf`s.
#' Therefore, `exact` is ignored in these cases and *p*-values of the respective
#' test settings are calculated by a normal approximation.
#'
#' By setting `exact = NULL`, exact computation is performed if the sample in a
#' test setting does not have any ties or zeros and if the sample size is lower
#' than or equal to 200.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`], [`signed_rank_test_pv()`],
#' [`mann_whitney_test_pv()`]
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
wilcox_single_test_pv <- function(
  x,
  mu = 0,
  alternative = "two.sided",
  exact = NULL,
  correct = TRUE,
  digits_rank = Inf,
  simple_output = FALSE
) {
  # plausibility checks of input parameters
  qassert(x, c("N+", "L+"))
  if(!is.list(x)) x <- list(x) else qassertr(x, "N+")
  len_x <- length(x)

  qassert(mu, "N+()")
  len_m <- length(mu)

  qassert(exact, c("B1", "0"))
  qassert(correct, "B1")

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
  W <- numeric(len_g)
  means <- numeric(len_g)
  sds <- numeric(len_g)
  zeros <- logical(len_g)
  ties <- logical(len_g)
  for(i in seq_len(len_g)) {
    y <- x[[i]] - mu[i]

    is_zero <- (y == 0)
    zeros[i] <- any(is_zero)
    if(zeros[i]) y <- y[!is_zero]

    n[i] <- length(y)

    ranks <- if(is.finite(digits_rank))
      rank(abs(signif(y, digits_rank))) else
        rank(abs(y))

    W[i] <- sum(ranks[y > 0])
    ties[i] <- length(ranks) != length(unique(ranks))

    means[i] <- n[i] * (n[i] + 1) / 4
    t <- table(ranks)
    sds[i] <- sqrt(n[i] * (n[i] + 1) * (2 * n[i] + 1) / 24 - sum(t^3 - t) / 48)
  }
  ex <- if(is.null(exact))
    !zeros & !ties & n < 201 else exact & !zeros & !ties & n < 1039

  # determine unique parameter sets
  params <- data.frame(alternative, n, ex, means, sds)
  params_ex <- unique(subset(params, ex, 1:2))
  params_ap <- unique(subset(params, !ex, -(2:3)))
  idx_ex   <- as.numeric(rownames(params_ex))
  idx_ap   <- as.numeric(rownames(params_ap))
  rows     <- c(idx_ex, idx_ap)
  params_u <- params[rows, -3]

  len_ex <- length(idx_ex)
  len_ap <- length(idx_ap)
  idx_ex <- seq_len(len_ex)
  idx_ap <- len_ex + seq_len(len_ap)
  len_u  <- len_ex + len_ap

  alts_u <- params_u$alternative
  n_u    <- params_u$n
  mean_u <- params_u$means
  sd_u   <- params_u$sds

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  if(!is.null(exact) && exact) {
    if(any(ties))
      warning("One or more p-values cannot be computed exactly because of ties")
    if(any(zeros))
      warning("One or more p-values cannot be computed exactly because of zeros")
    if(!any(ties) & !any(zeros) & any(n > 1038))
      warning(paste(
        "One or more p-values cannot be computed",
        "exactly because sample size exceeds 1,038"
      ))
  }

  # begin exact computations (if any)
  for(i in idx_ex) {
    idx_supp <- which(alts_u[i] == alternative & n_u[i] == n & ex)

    if(simple_output) {
      # compute p-values directly
      res[idx_supp] <- switch(
        EXPR = alts_u[i],
        less = psignrank(W[idx_supp], n_u[i]),
        greater = psignrank(W[idx_supp] - 1, n_u[i], lower.tail = FALSE),
        two.sided = {
          idx_l <- which(W[idx_supp] < mean_u[i])
          idx_u <- which(W[idx_supp] >= mean_u[i])
          pv <- numeric(length(idx_supp))
          if(length(idx_l))
            pv[idx_l] <- psignrank(W[idx_supp][idx_l], n_u[i])
          if(length(idx_u))
            pv[idx_u] <- psignrank(n_u[i] - W[idx_supp][idx_u], n_u[i])
          pmin(1, 2 * pv)
        }
      )
    } else {
      # generate all probabilities under current sample size
      d <- generate_signrank_probs(n_u[i])
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[i],
        probs = d,
        expectations = abs(seq_along(d) - 1 - mean_u[i])
      )

      # store results and support
      res[idx_supp] <- pv_supp[W[idx_supp] + 1]
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  # begin approximation computations (if any)
  for(i in idx_ap) {
    idx_supp <- which(alts_u[i] == alternative & !ex & mean_u[i] == means &
                        sd_u[i] == sds)

    if(simple_output) {
      z <- (W[idx_supp] - mean_u[i]) / sd_u[i]
      res[idx_supp] <- switch(
        EXPR = alts_u[i],
        less = ifelse(W[idx_supp] == 0, 0, pnorm(z + correct * 0.5 / sd_u[i])),
        greater = ifelse(
          W[idx_supp] == 0,
          1,
          pnorm(z - correct * 0.5 / sd_u[i], lower.tail = FALSE)
        ),
        two.sided = 2 * pnorm(-abs(z) + correct * 0.5 / sd_u[i])
      )
    } else {
      # compute p-value support
      pv_supp <- support_normal(
        alternative = alts_u[i],
        x = 0L:((n_u[i] * (n_u[i] + 1L)) %/% 2L),
        mean = mean_u[i],
        sd = sd_u[i],
        correct = correct
      )

      # store results and support
      res[idx_supp] <- pv_supp[W[idx_supp] + 1]
      if(!simple_output) {
        supports[[i]] <- unique(sort(pv_supp))
        indices[[i]]  <- idx_supp
      }
    }
  }

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon signed rank test",
      inputs = list(
        observations = x,
        nullvalues = data.frame(location = mu),
        parameters = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            n = ifelse(ex, n, NA),
            mean = ifelse(!ex, means, NA),
            sd = ifelse(!ex, sds, NA),
            `continuity correction` = ifelse(!ex, correct, NA),
            alternative = alternative,
            exact = ex,
            distribution = ifelse(ex, "Wilcoxon's signed-rank", "normal"),
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(W),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = dnames["x"]
    )
  } else res

  return(out)
}
