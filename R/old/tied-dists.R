d_MW_exact_ties <- function(scores, size_x) {
  scores <- sort(scores)

  min_scores <- cumsum(c(0L, scores[seq_len(size_x)]))
  max_scores <- cumsum(c(0L, rev(scores)[seq_len(size_x)]))

  size_l <- min(size_x, length(scores) - size_x)

  dp <- vector("list", size_l + 1)
  dp <- lapply(max_scores - min_scores + 1L, function(len) numeric(len))
  dp[[1]][1] <- 1
  first <- c(1, numeric(size_l))
  last  <- c(1, numeric(size_l))

  for(i in seq_along(scores)) {
    s <- scores[i]
    for(k in seq.int(min(i, size_l), 1, by = -1)) {
      w_keep <- (i - k) / i
      w_add  <- k / i

      if(first[k + 1]) {
        lo <- first[k + 1] - min_scores[k + 1]
        hi <- last[k + 1]  - min_scores[k + 1]
        dp[[k + 1]][lo:hi] <- dp[[k + 1]][lo:hi] * w_keep
      }

      src_lo <- first[k] - min_scores[k]
      src_hi <-  last[k] - min_scores[k]

      dst_lo <- first[k] - min_scores[k + 1] + s
      dst_hi <- last[k]  - min_scores[k + 1] + s

      dp[[k + 1]][dst_lo:dst_hi] <- dp[[k + 1]][dst_lo:dst_hi] +
        dp[[k]][src_lo:src_hi] * w_add

      new_min <- first[k] + s
      new_max <-  last[k] + s
      if(!first[k + 1]) {
        first[k + 1] <- new_min
        last[k + 1] <- new_max
      } else {
        if(new_min < first[k + 1]) first[k + 1] <- new_min
        if(new_max >  last[k + 1])  last[k + 1] <- new_max
      }
    }
  }

  counts <- dp[[size_l + 1]]
  probs  <- numerical_adjust(counts)

  if(size_x == size_l) return(probs) else return(rev(probs))
}

d_SR_exact_ties <- function(scores) {
  scores <- sort(scores)
  total  <- sum(scores)

  dp    <- numeric(total + 1)
  dp[1] <- 1
  last  <- 1L

  for (i in seq_along(scores)) {
    s <- scores[i]

    # alte Verteilung sichern, bevor sie überschrieben wird
    old <- dp[seq_len(last)]

    # Vorzeichen "-": alte Verteilung bleibt an derselben Stelle, Gewicht 1/2
    dp[seq_len(last)] <- dp[seq_len(last)] * 0.5

    # Vorzeichen "+": alte Verteilung wird um s verschoben, Gewicht 1/2
    dp[s + seq_len(last)] <- dp[s + seq_len(last)] + old * 0.5

    last <- last + s
  }

  probs <- dp[seq_len(last)]
  probs <- numerical_adjust(probs)

  return(probs)
}
