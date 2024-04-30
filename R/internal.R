#'@importFrom stats stepfun
ts.pv <- function(statistics, probs, decreasing = FALSE, normalize = FALSE){
  ord.stats <- order(statistics, decreasing = decreasing)
  statistics <- statistics[ord.stats]
  probs <- probs[ord.stats]
  pk <- cumsum(probs)
  pk[length(pk)] <- 1
  ord.sf <- if(decreasing) order(statistics) else 1:length(statistics)
  sf <- stepfun(statistics[ord.sf], c(0, pk[ord.sf]))
  pv <- sf(statistics)

  if(normalize) pv <- pv/max(pv)
  return(pv[order(ord.stats)])
}

numerical.adjust <- function(P, normalize = TRUE, rel.tol = .Machine$double.eps * 128){
  ord <- order(P)
  P <- P[ord]
  idx.dup <- which(duplicated(P))
  Q <- unique(P)
  n <- length(Q) - 1
  rel.diff <- abs(Q[1:n] / Q[1:n + 1] - 1)
  idx.diff <- which(rel.diff <= rel.tol & rel.diff > 0)
  len <- length(idx.diff)
  if(len){
    for(i in 1:len){
      j <- 1
      while(i < len && i + j <= len && idx.diff[i] + j <= n && idx.diff[i + j] == idx.diff[i] + j) j <- j + 1
      Q[idx.diff[i]:(idx.diff[i] + j)] <- mean(Q[idx.diff[i]:(idx.diff[i] + j)])
    }
  }
  if(length(idx.dup)){
    P[-idx.dup] <- Q
    P[idx.dup] <- -Inf
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
generate.binom.probs <- function(n, p, log = FALSE){
  d <- numerical.adjust(dbinom(0:n, n, p))

  if(log) return(log(d)) else return(d)
}

#'@importFrom stats dpois qpois
generate.poisson.probs <- function(lambda, log = FALSE){
  # search for last observation with P(X = limit) > 0
  limit <- ifelse(lambda == 0, 0, qpois(2^-1074, lambda, FALSE))
  d <- dpois(0:limit, lambda, log)

  return(numerical.adjust(d))
}

support_exact <- function(alternative, probs, expectations = NULL){
  support <- pmin(
    1,
    switch(
      EXPR    = alternative,
      less    = c(cumsum(probs[-length(probs)]), 1),
      greater = c(1, rev(cumsum(rev(probs[-1])))),
      minlike = ts.pv(statistics = probs, probs = probs),
      blaker  = ts.pv(
        statistics = pmin(
          c(cumsum(probs[-length(probs)]), 1),
          c(1, rev(cumsum(rev(probs[-1]))))
        ),
        probs = probs
      ),
      absdist = ts.pv(
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
pnorm_zero <- function(q, sd = 1, lower.tail = TRUE){
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
    if(lower.tail)
      res[idx0] <- as.numeric(q[idx0] >= 0) else
        res[idx0] <- as.numeric(q[idx0] <= 0)
  }
  if(length(idx1))
    res[idx1] <- pnorm(q[idx1], sd = sd[idx1], lower.tail = lower.tail)

  # return results
  return(res)
}
