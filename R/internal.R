#' @importFrom stats stepfun
pvalues_by_statistics <- function(statistics, probs, decreasing = FALSE, normalize = FALSE){
  stats_order <- order(statistics, decreasing = decreasing)
  statistics <- statistics[stats_order]
  probs <- probs[stats_order]
  probs_cumulated <- cumsum(probs)
  probs_cumulated[length(probs_cumulated)] <- 1
  order_for_stepfun <- if(decreasing) order(statistics) else seq_along(statistics)
  probs_stepfun <- stepfun(statistics[order_for_stepfun],
                           c(0, probs_cumulated[order_for_stepfun]))
  pvalues <- probs_stepfun(statistics)

  if(normalize) pvalues <- pvalues/max(pvalues)
  return(pvalues[order(stats_order)])
}

numerical_adjust <- function(values, normalize = TRUE, rel.tol = .Machine$double.eps * 128) {
  values_order <- order(values)
  values <- values[values_order]
  idx_duplicates <- which(duplicated(values))
  values_unique <- unique(values)
  len_values_unique <- length(values_unique) - 1
  rel_diffs <- (values_unique[1:len_values_unique + 1] - values_unique[1:len_values_unique]) / values_unique[1:len_values_unique + 1]
  idx_diffs <- which(rel_diffs <= rel.tol & rel_diffs > 0)
  len_diffs <- length(idx_diffs)
  for(i in seq_len(len_diffs)) {
    j <- 1

    while(
      i < len_diffs && i + j <= len_diffs &&
      idx_diffs[i] + j <= len_values_unique &&
      idx_diffs[i + j] == idx_diffs[i] + j
    ) j <- j + 1

    values_unique[idx_diffs[i] + 0:j] <- mean(values_unique[idx_diffs[i] + 0:j])
  }
  if(length(idx_duplicates)) {
    values[-idx_duplicates] <- values_unique
    values[idx_duplicates] <- -Inf
    values <- cummax(values)
  } else values <- values_unique
  values <- values[order(values_order)]

  if(normalize) {
    sum_values <- sum(values)
    if(sum_values <= 0){
      values <- exp(values)
      sum_values <- sum(values)
    }
    counter <- 1
    while(sum_values != 1 && (sum_values > 1 || counter <= 10)) {
      values <- values/sum_values
      sum_values <- sum(values)
      counter <- counter + 1
    }
    #if(sum_values > 1 && counter > 10)
    #  warning("Sum of probabilities slightly exceeds 1. Normalisation attempt failed.\n")
  }

  return(values)
}

#' @importFrom stats dbinom
generate_binom_probs <- function(size, prob, log = FALSE){
  probability_masses <- numerical_adjust(dbinom(0:size, size, prob))

  if(log) return(log(probability_masses)) else return(probability_masses)
}

#' @importFrom stats dpois qpois
generate_poisson_probs <- function(lambda, log = FALSE){
  # search for last observation with P(X = limit) > 0
  limit <- ifelse(lambda == 0, 0, qpois(2^-1074, lambda, FALSE))
  probability_masses <- numerical_adjust(dpois(0:limit, lambda))

  if(log) return(log(probability_masses)) else return(probability_masses)
}

#' @importFrom stats dsignrank
generate_sign_rank_probs <- function(n, log = FALSE) {
  len <- length(n)
  if(len == 0) return(list())

  ord_n <- order(n)
  n <- n[ord_n]
  max_n <- n[len]

  res <- sign_rank_probs_int(n)[order(ord_n)]

  for(i in seq_len(len)) {
    res[[i]] <- if(log)
      log(numerical_adjust(res[[i]])) else
        numerical_adjust(res[[i]])
  }

  return(res)
}

generate_mann_whitney_probs <- function(m, n, log = FALSE) {
  len <- length(m)
  if(len == 0) return(list())
  sizes <- cbind(pmax(m, n), pmin(m, n))

  ord_n <- order(sizes[, 2])
  sizes <- sizes[ord_n, , drop = FALSE]
  ord_m <- order(sizes[, 1])
  sizes <- sizes[ord_m, , drop = FALSE]
  m <- as.integer(sizes[, 1])
  n <- as.integer(sizes[, 2])

  res <- mann_whitney_probs_int(m, n)[order(ord_m)][order(ord_n)]

  for(i in seq_len(len)) {
    res[[i]] <- if(log)
      log(numerical_adjust(res[[i]])) else
        numerical_adjust(res[[i]])
  }

  return(res)
}

support_exact <- function(alternative, probs, expectations = NULL){
  support <- pmin(
    1,
    switch(
      EXPR    = alternative,
      less    = c(cumsum(probs[-length(probs)]), 1),
      greater = c(1, rev(cumsum(rev(probs[-1])))),
      minlike = pvalues_by_statistics(statistics = probs, probs = probs),
      blaker  = pvalues_by_statistics(
        statistics = pmin(cumsum(probs), rev(cumsum(rev(probs)))),
        probs = probs
      ),
      absdist = pvalues_by_statistics(
        statistics = expectations,
        probs = probs,
        decreasing = TRUE
      ),
      # default for "two.sided" or "central", i.e. double of smaller tail
      2 * pmin(cumsum(probs), rev(cumsum(rev(probs))))
    )
  )
  return(support)
}

#' @importFrom stats pnorm
support_normal <- function(alternative, x, mean, sd, correct) {
  z <- (x - mean) / sd
  support <- pmin(
    1,
    switch(
      EXPR      = alternative,
      less      = rev(c(1, pnorm(rev(z)[-1] + correct * 0.5 / sd))),
      greater   = c(1, pnorm(z[-1] - correct * 0.5 / sd, lower.tail = FALSE)),
      # defaulting for "two.sided"
      2 * pnorm(-abs(z) + correct * 0.5 / sd)
    )
  )

  return(support)
}

# modification of pnorm that ensures that P(X >= q) is computed for sd = 0 (instead of P(X > q))
#' @importFrom stats pnorm
pnorm_zero <- function(q, sd = 1, lower_tail = TRUE){
  # number of input values
  len <- length(q)
  # vector of results
  res <- numeric(len)
  # make sure 'sd' has the same size as 'q'
  sd <- rep_len(sd, len)
  # 'sd = 0' means we have a single-point distribution at 0
  idx0 <- which(sd == 0)
  # standard normal distribution
  idx1 <- which(sd != 0)
  if(length(idx0)){
    if(lower_tail)
      res[idx0] <- as.numeric(q[idx0] >= 0) else
        res[idx0] <- as.numeric(q[idx0] <= 0)
  }
  if(length(idx1))
    res[idx1] <- pnorm(q[idx1], sd = sd[idx1], lower.tail = lower_tail)

  # return results
  return(res)
}

p_from_d <- function(q, pmf, lower.tail = TRUE) {
  res <- if(lower.tail) {
    cumsum(pmf[seq_len(max(q) + 1)])[q + 1]
  } else {
    rev(cumsum(rev(pmf[(min(q) + 2):length(pmf)])))[q - min(q) + 1]
  }
  return(res)
}
