#' @name mann_whitney_test_pv
#'
#' @title
#' Wilcoxon-Mann-Whitney *U* test
#'
#' @description
#' `mann_whitney_test_pv()` performs an exact or approximate
#' Wilcoxon-Mann-Whitney *U* test about the location shift between two
#' independent groups when the data is not necessarily normally distributed. In
#' contrast to [`stats::wilcox.test()`], it is vectorised and only calculates
#' *p*-values. Furthermore, it is capable of returning the discrete *p*-value
#' supports, i.e. all observable *p*-values under a null hypothesis. Multiple
#' tests can be evaluated simultaneously.
#'
#' @param x,y     numerical vectors forming the samples to be tested or lists
#'                of numerical vectors for multiple tests.
#' @param mu      numerical vector of hypothesised location shift(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' We use a test statistic called the Wilcoxon Rank Sum Statistic, defined by
#' \deqn{U = \sum_{i = 1}^{n_X}{rank(X_i)} - \frac{n_X(n_X + 1)}{2},}
#' where \eqn{rank(X_i)} is the rank of \eqn{X_i} in the concatenated sample
#' of \eqn{X} and \eqn{Y}, and \eqn{n_X} and \eqn{n_Y} are the respective
#' sizes of the samples \eqn{X} and \eqn{Y}. Note that \eqn{U}
#' can range from \eqn{0} to \eqn{n_X \cdot n_Y}.
#' This is the same statistic used by [`stats::wilcox.test()`] and
#' whose distribution is accessible with [`pwilcox`].
#' This is also the statistic defined by the two given references.
#' Note, however, that it is not what is called the Mann-Whitney U Statistic
#' in the (English-language) Wikipedia article (as of February 12, 2026). The
#' latter is defined as, using our notation, \eqn{\min(U, n_X \cdot n_Y - U)}.
#' Using the Wikipedia notation, the Wilcoxon Rank Sum Statistic is \eqn{U_2}.
#'
#' The parameters `x`, `y`, `mu` and `alternative` are vectorised. If `x` and
#' `y` are lists, they are replicated automatically to have the same lengths. In
#' case `x` or `y` are not lists, they are added to new ones, which are then
#' replicated to the appropriate lengths. This allows multiple hypotheses to be
#' tested simultaneously.
#'
#' In the presence of ties, computation of exact *p*-values is not possible.
#' Therefore, `exact` is ignored in these cases and *p*-values of the
#' respective test settings are calculated by a normal approximation.
#'
#' By setting `exact = NULL`, exact computation is performed only if both
#' samples in a test setting do not have any ties and if both sample sizes are
#' lower than or equal to 200. If any of these conditions is not met,
#' \eqn{p}-values are computed by normal approximation.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`], [`pwilcox`], [`wilcox_test_pv()`]
#'
#' @references
#' Mann, H. D. & Whitney, D. R. (1947). On a Test of Whether one of Two Random
#'   Variables is Stochastically Larger than the Other. *Ann. Math. Statist.*,
#'   *18*(1), pp. 50-60. \doi{10.1214/aoms/1177730491}
#'
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. pp. 115-135. \doi{10.1002/9781119196037}
#'
#' @examples
#' # Constructing
#' set.seed(1)
#' r1 <- rnorm(100)
#' r2 <- rnorm(100, 1)
#'
#' # Exact two-sided p-values and their supports
#' results_ex  <- mann_whitney_test_pv(r1, r2)
#' print(results_ex)
#' results_ex$get_pvalues()
#' results_ex$get_pvalue_supports()
#'
#' # Normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- mann_whitney_test_pv(r1, r2, alternative = "less", exact = FALSE)
#' print(results_ap)
#' results_ap$get_pvalues()
#' results_ap$get_pvalue_supports()
#'
#' @importFrom checkmate qassert qassertr
#' @importFrom cli cli_warn
#' @export
mann_whitney_test_pv <- function(
  x,
  y,
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

  qassert(y, c("N+", "L+"))
  if(!is.list(y)) y <- list(y) else qassertr(y, "N+")
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

  qassert(digits_rank, "N")

  qassert(simple_output, "B1")

  # replicate inputs to same length
  len_g <- max(len_x, len_y, len_m, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_y < len_g) y <- rep_len(y, len_g)
  if(len_m < len_g) mu <- rep_len(mu, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # compute ranks and lengths
  nx    <- integer(len_g)
  ny    <- integer(len_g)
  N     <- integer(len_g)
  nn    <- integer(len_g)
  U     <- numeric(len_g)
  means <- numeric(len_g)
  sds   <- numeric(len_g)
  ties  <- logical(len_g)
  stats <- vector("list", len_g)
  for(i in seq_len(len_g)) {
    nx[i] <- length(x[[i]])
    ny[i] <- length(y[[i]])
    N[i]  <- nx[i] + ny[i]
    nn[i] <- nx[i] * ny[i]

    ranks <- if(is.finite(digits_rank))
      rank(signif(c(x[[i]] - mu[i], y[[i]]), digits_rank)) else
        rank(c(x[[i]] - mu[i], y[[i]]))

    U[i] <- sum(ranks[seq_len(nx[i])]) - nx[i] * (nx[i] + 1) / 2
    ties[i] <- length(ranks) != length(unique(ranks))
    if(ties[i]) stats[[i]] <- as.integer(round(2 * ranks))

    means[i] <- nn[i] / 2
    t <- table(ranks)
    sds[i] <- sqrt((nn[i] / 12) * ((N[i] + 1) -
                      sum(t^3 - t)/(N[i] * (N[i] - 1))))
  }

  ex <- if(is.null(exact)) nx < 201 & ny < 201 else rep(exact, len_g)
  ew <- edgeworth & !ex & !ties
  ew[ex | ties] <- NA

  # compute Edgeworth coefficients for normal approximations, if desired
  idx_ew <- which(ew)
  if(length(idx_ew)) {
    ew_coefs <- matrix(NA_real_, len_g, edgeworth)
    if(edgeworth >= 1)
      ew_coefs[idx_ew, 1] <- -(N[idx_ew]^2 - nn[idx_ew] + N[idx_ew]) /
        (20 * nn[idx_ew] * (N[idx_ew] + 1))
    if(edgeworth >= 2)
      ew_coefs[idx_ew, 2] <- (
        2*(nx[idx_ew]^4 + ny[idx_ew]^4) +
        4*(nn[i]*(nx[idx_ew]^2 + ny[idx_ew]^2) + nx[idx_ew]^3 + ny[idx_ew]^3) +
        6*nx[idx_ew]^2 * ny[idx_ew]^2 +
        N[idx_ew] * (7*nn[idx_ew] + N[idx_ew] - 1)
      ) / (210 * nx[idx_ew]^2 * ny[idx_ew]^2 * (N[idx_ew] + 1)^2)
    if(edgeworth == 3)
      ew_coefs[idx_ew, 3] <- ew_coefs[idx_ew, 1]^2 / 2
  }

  # determine unique parameter sets
  params    <- data.frame(alternative, nx, ny, means, sds, ew, ties)
  params_ex <- unique(subset(params, ex & !ties, 1:3))
  params_ti <- subset(params, ex & ties)
  params_ap <- unique(subset(params, !ex, -(2:3)))
  idx_ex    <- as.numeric(rownames(params_ex))
  idx_ti    <- as.numeric(rownames(params_ti))
  idx_ap    <- as.numeric(rownames(params_ap))
  rows      <- c(idx_ex, idx_ti, idx_ap)
  params_u  <- params[rows, ]
  if(any(ew, na.rm = TRUE)) ew_coefs_u <- ew_coefs[rows, , drop = FALSE]

  len_ex <- length(idx_ex)
  len_ti <- length(idx_ti)
  len_ap <- length(idx_ap)
  idx_ex <- seq_len(len_ex)
  idx_ti <- len_ex + seq_len(len_ti)
  idx_ap <- len_ex + len_ti + seq_len(len_ap)
  len_u  <- len_ex + len_ti + len_ap

  alts_u <- params_u$alternative
  nx_u   <- params_u$nx
  ny_u   <- params_u$ny
  mean_u <- params_u$means
  sd_u   <- params_u$sds
  ew_u   <- params_u$ew
  ties_u <- params_u$ties

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # pre-compute exact distributions (if any)
  sizes_ex <- unique(data.frame(nx_u, ny_u)[idx_ex, ])
  d <- generate_mann_whitney_probs(sizes_ex[, 1], sizes_ex[, 2])

  # begin exact computations for settings without ties (if any)
  for(i in idx_ex) {
    idx_par <- which(
      alts_u[i] == alternative & nx_u[i] == nx & ny_u[i] == ny & ex & !ties
    )

    idx_d <- which(sizes_ex[, 1] == nx_u[i] & sizes_ex[, 2] == ny_u[i])

    if(simple_output) {
      # compute p-values directly
      res[idx_par] <- switch(
        EXPR = alts_u[i],
        less = p_from_d(U[idx_par], d[[idx_d]]),
        greater = p_from_d(U[idx_par] - 1, d[[idx_d]], FALSE),
        two.sided = {
          idx_l <- which(U[idx_par] < mean_u[i])
          idx_u <- which(U[idx_par] >= mean_u[i])
          pv <- numeric(length(idx_par))
          if(length(idx_l))
            pv[idx_l] <- p_from_d(U[idx_par][idx_l], d[[idx_d]])
          if(length(idx_u))
            pv[idx_u] <- p_from_d(U[idx_par][idx_u] - 1, d[[idx_d]], FALSE)
          pmin(1, 2 * pv)
        }
      )
    } else {
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[i],
        probs = d[[idx_d]]
      )

      # store results and support
      res[idx_par] <- pv_supp[U[idx_par] + 1]
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_par
    }
  }

  # begin exact computations for settings with ties (if any)
  idx_par <- which(ex & ties)
  S       <- U[idx_par]
  stats   <- stats[idx_par]
  for(i in seq_along(idx_par)) {
    j <- idx_ti[i]

    # compute exact distribution with ties
    dist <- numerical_adjust(mann_whitney_probs_ties_int(stats[[i]], nx_u[j]))
    min_stat <- sum(sort(stats[[i]])[seq_len(nx_u[j])])
    R <- 2 * S[i] + nx_u[j] * (nx_u[j] + 1) - min_stat

    if(simple_output) {
      # compute p-values directly
      res[idx_par[i]] <- switch(
        EXPR = alts_u[j],
        less = p_from_d(R, dist),
        greater = p_from_d(R - 1, dist, FALSE),
        two.sided = pmin(1, 2 * p_from_d(
          R - (S[i] > mean_u[j]), dist, S[i] <= mean_u[j]
        ))
      )
    } else {
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[j],
        probs = dist
      )

      # store results and support
      res[idx_par[i]] <- pv_supp[R + 1]
      supports[[j]] <- unique(sort(pv_supp))
      indices[[j]]  <- idx_par[i]
    }
  }

  # begin approximation computations (if any)
  for(i in idx_ap) {
    idx_par <- which(
      alts_u[i] == alternative & !ex & ties_u[i] == ties &
        mean_u[i] == means & sd_u[i] == sds
    )

    ew_ok <- !is.na(ew_u[i]) && ew_u[i] && !ties_u[i]
    e <- if(ew_ok) ew_coefs_u[i, ]

    if(simple_output) {
      res[idx_par] <- switch(
        EXPR = alts_u[i],
        less = pnorm_MW_edgeworth(
          U[idx_par], mean_u[i], sd_u[i], TRUE, correct, e
        ),
        greater = pnorm_MW_edgeworth(
          U[idx_par], mean_u[i], sd_u[i], FALSE, correct, e
        ),
        two.sided = pmin(1, 2 * if(ew_ok)
          pmin(
            pnorm_MW_edgeworth(
              U[idx_par], mean_u[i], sd_u[i], TRUE, correct, e
            ),
            pnorm_MW_edgeworth(
              U[idx_par], mean_u[i], sd_u[i], FALSE, correct, e
            )
          ) else pmin(1,
            pnorm(-abs(U[idx_par] - mean_u[i]), -correct * 0.5, sd_u[i])
          )
        )
      )
    } else {
      # compute p-value support
      z <- if(!any(ties[idx_par])) 0L:(nx_u[i] * ny_u[i]) else
        seq(0, nx_u[i] * ny_u[i], 0.5)
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
      idx_supp <- 1 + if(!any(ties[idx_par]))
        U[idx_par] else round(2 * U[idx_par])
      res[idx_par] <- pv_supp[idx_supp]
      if(!simple_output) {
        supports[[i]] <- unique(sort(pv_supp))
        indices[[i]]  <- idx_par
      }
    }
  }

  # create output object
  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon-Mann-Whitney U test",
      inputs = list(
        observations = list(x, y),
        parameters = NULL,
        nullvalues = data.frame(`location shift` = mu, check.names = FALSE),
        computation = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            alternative = alternative,
            exact = ex,
            distribution = ifelse(
              ex,
              paste0(
                "Wilcoxon-Mann-Whitney", ifelse(ties, " (tie-adjusted)", "")
              ),
              "normal"
            ),
            #distribution.mean = ifelse(!ex, means, NA_real_),
            #distribution.sd = ifelse(!ex, sds, NA_real_),
            `continuity correction` = ifelse(ex, NA, correct),
            `Edgeworth expansion` = ew,
            `Edgeworth series terms` = ifelse(ew, ew * edgeworth, NA),
            `size of first sample` = nx,
            `size of second sample` = ny,
            ties = ties,
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(U),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = dnames[c("x", "y")]
    )
  } else res

  # return results
  return(out)
}
