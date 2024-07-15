#'@importFrom stats stepfun
ts_pv <- function(statistics, probs, decreasing = FALSE, normalize = FALSE){
  ord_stats <- order(statistics, decreasing = decreasing)
  statistics <- statistics[ord_stats]
  probs <- probs[ord_stats]
  pk <- cumsum(probs)
  pk[length(pk)] <- 1
  ord_sf <- if(decreasing) order(statistics) else 1:length(statistics)
  sf <- stepfun(statistics[ord_sf], c(0, pk[ord_sf]))
  pv <- sf(statistics)

  if(normalize) pv <- pv/max(pv)
  return(pv[order(ord_stats)])
}

numerical_adjust <- function(P, normalize = TRUE, rel.tol = .Machine$double.eps * 128){
  ord <- order(P)
  P <- P[ord]
  idx_dup <- which(duplicated(P))
  Q <- unique(P)
  n <- length(Q) - 1
  rel_diff <- abs(Q[1:n] / Q[1:n + 1] - 1)
  idx_diff <- which(rel_diff <= rel.tol & rel_diff > 0)
  len <- length(idx_diff)
  if(len){
    for(i in 1:len){
      j <- 1
      while(i < len && i + j <= len && idx_diff[i] + j <= n && idx_diff[i + j] == idx_diff[i] + j) j <- j + 1
      Q[idx_diff[i]:(idx_diff[i] + j)] <- mean(Q[idx_diff[i]:(idx_diff[i] + j)])
    }
  }
  if(length(idx_dup)){
    P[-idx_dup] <- Q
    P[idx_dup] <- -Inf
    P <- cummax(P)
  } else P <- Q
  P <- P[order(ord)]

  if(normalize){
    s <- sum(P)
    if(s <= 0){
      P <- exp(P)
      s <- sum(P)
    }
    k <- 1
    while(s != 1 && (s > 1 || k <= 10)){
      P <- P/s
      s <- sum(P)
      k <- k + 1
    }
    if(s > 1 && k > 10)
      warning("Sum of probabilities slightly exceeds 1. Normalisation attempt failed.\n")
  }

  return(P)
}

#'@importFrom stats dbinom
generate_binom_probs <- function(n, p, log = FALSE){
  d <- numerical_adjust(dbinom(0:n, n, p))

  if(log) return(log(d)) else return(d)
}

#'@importFrom stats dpois qpois
generate_poisson_probs <- function(lambda, log = FALSE){
  # search for last observation with P(X = limit) > 0
  limit <- ifelse(lambda == 0, 0, qpois(2^-1074, lambda, FALSE))
  d <- dpois(0:limit, lambda, log)

  return(numerical_adjust(d))
}

support_exact <- function(alternative, probs, expectations = NULL){
  support <- pmin(
    1,
    switch(
      EXPR    = alternative,
      less    = c(cumsum(probs[-length(probs)]), 1),
      greater = c(1, rev(cumsum(rev(probs[-1])))),
      minlike = ts_pv(statistics = probs, probs = probs),
      blaker  = ts_pv(
        statistics = pmin(
          c(cumsum(probs[-length(probs)]), 1),
          c(1, rev(cumsum(rev(probs[-1]))))
        ),
        probs = probs
      ),
      absdist = ts_pv(
        statistics = expectations,
        probs = probs,
        decreasing = TRUE
      ),
      central = 2 * pmin(
        c(cumsum(probs[-length(probs)]), 1),
        c(1, rev(cumsum(rev(probs[-1]))))
      )
    )
  )
  return(support)
}

# modification of pnorm that ensures that P(X >= q) is computed for sd = 0 (instead of P(X > q))
#'@importFrom stats pnorm
pnorm_zero <- function(q, sd = 1, lower_tail = TRUE){
  # number of input values
  n <- length(q)
  # vector of results
  res <- numeric(n)
  # make sure 'sd' has the same size as 'q'
  sd <- rep_len(sd, n)
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
