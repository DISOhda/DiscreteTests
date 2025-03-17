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
#' @param x,y           numerical vectors forming the samples to be tested or
#'                      lists of numerical vectors for multiple samples.
#' @param d             numerical vector of hypothesised differences(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#' @templateVar digits_rank TRUE
#'
#' @details
#' The parameters `x`, `y`, `d` and `alternative` are vectorised. If `x` and `y`
#' are lists, they are replicated automatically to have the same lengths. In
#' case `x` or `y` are not lists, they are added to new ones, which are then
#' replicated to the appropriate lengths. This allows multiple hypotheses to be
#' tested simultaneously.
#'
#' In the presence of ties, computation of exact *p*-values is not possible.
#' Therefore, `exact` is ignored in these cases and *p*-values of the
#' respective test settings are calculated by a normal approximation.
#'
#' Similarly, if the sum of the sample sizes of `x` and `y` is greater than
#' 1,000, [`stats::dwilcox()`] tends to produce `NaN`s or only zeros.
#' Additionally, the memory requirements may exceed the available RAM of the
#' user's system. The user should therefore avoid exact computation with large
#' sample sizes, as it provides only a minor gain in accuracy over the normal
#' approximation.
#'
#' By setting `exact = NULL`, exact computation is performed if both samples in
#' a test setting do not have any ties and if *both* sample sizes are lower than
#' or equal to 200.
#'
#' If `digits_rank = Inf` (the default), [`rank()`][`base::rank()`] is used to
#' compute ranks for the tests statistics instead of
#' [`rank`][`base::rank()`]([`signif(., digits_rank)`][`base::signif()`])
#'
#' @template return
#'
#' @seealso
#' [`stats::wilcox.test()`], [`wilcox_single_test_pv()`],
#' [`signed_rank_test_pv()`]
#'
#' @references
#' Hollander, M. & Wolfe, D. (1973). *Nonparametric Statistical Methods*. Third
#'   Edition. New York: Wiley. pp. 115-125. \doi{10.1002/9781119196037}
#'
#' @examples
#' # Constructing
#' set.seed(1)
#' r1 <- rnorm(100)
#' r2 <- rnorm(100, 1)
#'
#' # Computation of exact two-sided p-values and their supports
#' results_ex  <- mann_whitney_test_pv(r1, r2)
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- mann_whitney_test_pv(r1, r2, alternative = "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom stats pwilcox
#' @importFrom checkmate qassert qassertr
#' @export
mann_whitney_test_pv <- function(
  x,
  y,
  d = 0,
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

  qassert(d, "N+()")
  len_d <- length(d)

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
  len_g <- max(len_x, len_y, len_d, len_a)
  if(len_x < len_g) x <- rep_len(x, len_g)
  if(len_y < len_g) y <- rep_len(y, len_g)
  if(len_d < len_g) d <- rep_len(d, len_g)
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  # compute ranks and lengths
  nx <- integer(len_g)
  ny <- integer(len_g)
  U <- numeric(len_g)
  means <- numeric(len_g)
  vars <- numeric(len_g)
  ties <- logical(len_g)
  for(i in seq_len(len_g)) {
    nx[i] <- length(x[[i]])
    ny[i] <- length(y[[i]])

    ranks <- if(is.finite(digits_rank))
      rank(signif(c(x[[i]] - d[i], y[[i]]), digits_rank)) else
        rank(c(x[[i]] - d[i], y[[i]]))

    U[i] <- sum(ranks[seq_len(nx[i])]) - nx[i] * (nx[i] + 1) / 2
    ties[i] <- length(ranks) != length(unique(ranks))

    means[i] <- nx[i] * ny[i] / 2
    t <- table(ranks)
    vars[i] <- sqrt((nx[i] * ny[i] / 12) * ((nx[i] + ny[i] + 1) -
                       sum(t^3 - t)/((nx[i] + ny[i]) * (nx[i] + ny[i] - 1))))
  }
  ex <- if(is.null(exact)) !ties & nx < 201 & ny < 201 else exact & !ties

  # determine unique parameter sets
  params <- data.frame(alternative, nx, ny, ex, means, vars)
  params_ex <- unique(subset(params, ex, 1:3))
  params_ap <- unique(subset(params, !ex, -(2:4)))
  idx_ex   <- as.numeric(rownames(params_ex))
  idx_ap   <- as.numeric(rownames(params_ap))
  rows     <- c(idx_ex, idx_ap)
  params_u <- params[rows, ]

  len_ex <- length(idx_ex)
  len_ap <- length(idx_ap)
  idx_ex <- seq_len(len_ex)
  idx_ap <- len_ex + seq_len(len_ap)
  len_u  <- len_ex + len_ap

  alts_u  <- params_u$alternative
  nx_u    <- params_u$nx
  ny_u    <- params_u$ny
  means_u <- params_u$means
  vars_u  <- params_u$vars

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  if(!is.null(exact) && exact) {
    if(any(ties))
      warning("One or more p-values cannot be computed exactly because of ties")
    if(any(nx + ny > 1000))
      warning(paste(
        "One or more p-values cannot be computed exactly",
        "because sums of sample pair sizes exceed 1,000"
      ))
  }

  # begin exact computations (if any)
  for(i in idx_ex) {
    idx_par <- which(alts_u[i] == alternative & nx_u[i] == nx & ny_u[i] == ny &
                       ex)

    if(simple_output) {
      # compute p-values directly
      res[idx_par] <- switch(
        EXPR = alts_u[i],
        less = pwilcox(U[idx_par], nx_u[i], ny_u[i]),
        greater = pwilcox(U[idx_par] - 1, nx_u[i], ny_u[i], lower.tail = FALSE),
        two.sided = {
          idx_l <- which(U[idx_par] < means_u[i])
          idx_u <- which(U[idx_par] >= means_u[i])
          pv <- numeric(length(idx_par))
          if(length(idx_l))
            pv[idx_l] <- pwilcox(U[idx_par][idx_l], nx_u[i], ny_u[i])
          if(length(idx_u))
            pv[idx_u] <- pwilcox(nx_u[i] * ny_u[i] - U[idx_par][idx_u],
                                 nx_u[i], ny_u[i])
          pmin(1, 2 * pv)
        }
      )
    } else {
      # generate all probabilities under current sample sizes
      probs <- generate_wilcox_probs(nx_u[i], ny_u[i])
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alts_u[i],
        probs = probs,
        expectations = abs(seq_along(d) - 1 - means_u[i])
      )

      # store results and support
      res[idx_par] <- pv_supp[U[idx_par] + 1]
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_par
    }
  }

  # begin approximation computations (if any)
  for(i in idx_ap) {
      idx_par <- which(alts_u[i] == alternative & !ex & means_u[i] == means &
                         vars_u[i] == vars)

    if(simple_output) {
      res[idx_par] <- support_normal(
        alternative = alts_u[i],
        x = U[idx_par],
        mean = means_u[i],
        sd = vars_u[i],
        correct = correct
      )
    } else {
      # compute p-value support
      pv_supp <- support_normal(
        alternative = alts_u[i],
        x = 0L:(nx_u[i] * ny_u[i]),
        mean = means_u[i],
        sd = vars_u[i],
        correct = correct
      )

      # store results and support
      res[idx_par] <- pv_supp[U[idx_par] + 1]
      if(!simple_output) {
        supports[[i]] <- unique(sort(pv_supp))
        indices[[i]]  <- idx_par
      }
    }
  }

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)

    DiscreteTestResults$new(
      test_name = "Wilcoxon-Mann-Whitney U test",
      inputs = list(
        observations = list(x, y),
        nullvalues = data.frame(`location shift` = d, check.names = FALSE),
        parameters = Filter(
          function(df) !all(is.na(df)),
          data.frame(
            `first sample size` = ifelse(ex, nx, NA),
            `second sample size` = ifelse(ex, ny, NA),
            mean = ifelse(!ex, means, NA),
            sd = ifelse(!ex, sqrt(vars), NA),
            alternative = alternative,
            exact = ex,
            distribution = ifelse(ex, "Wilcoxon-Mann-Whitney", "normal"),
            check.names = FALSE
          )
        )
      ),
      statistics = data.frame(U),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = paste(dnames["x"], "and", dnames["y"])
    )
  } else res

  return(out)
}
