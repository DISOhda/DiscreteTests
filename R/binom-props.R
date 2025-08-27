#' @title Unconditional two-sample Binomial test
#'
#' @export
binom_prop_test_pv <- function(
  x,
  n,
  p = NULL,
  alternative = "two.sided",
  #method = "pooled",
  pooled = TRUE,
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # plausibility checks of input parameters

  # define error message for malformed x or n
  error_msg <- paste("must either be a two-element vector, a two-column matrix",
                     "or a two-column data frame")

  # make sure that x is a data frame, matrix or vector
  qassert(x, c("D+", "M+", "V2"))
  # make sure that n is a data frame, matrix or vector
  qassert(n, c("D+", "M+", "V2"))
  # make sure that p is a data frame, matrix, vector or NULL
  qassert(p, c("D+", "M+", "V<=2", "0"))

  # if x is an atomic vector, make it a matrix with one row
  if(is.vector(x) && is.atomic(x))
    x <- t(x)
  # if n is an atomic vector, make it a matrix with one row
  if(is.vector(n) && is.atomic(n))
    n <- t(n)
  # if p is an atomic vector, make it a matrix with one row
  if(is.vector(p) && is.atomic(p))
    p <- t(p)

  # if x is a data frame, make it a matrix
  if(is.data.frame(x))
    x <- as.matrix(x)
  # if n is a data frame, make it a matrix
  if(is.data.frame(n))
    n <- as.matrix(n)
  # if p is a data frame, make it a matrix
  if(is.data.frame(p))
    p <- as.matrix(p)

  # when x is a matrix, it must satisfy some conditions
  if(is.matrix(x)) {
    # stop immediately, if dimensions are violated
    if(ncol(x) != 2)
      cli_abort(paste("'x'", error_msg))
    # check if all values are non-negative and close to integer
    assert_integerish(x, lower = 0, min.len = 2)
    # coerce to integer
    x <- round(x)
    mode(x) <- "integer"
  } else cli_abort(paste("'x'", error_msg))
  # when n is a matrix, it must satisfy some conditions
  if(is.matrix(n)) {
    # stop immediately, if dimensions are violated
    if(ncol(n) != 2)
      cli_abort(paste("'n'", error_msg))
    # check if all values are non-negative and close to integer
    assert_integerish(n, lower = 0, min.len = 2)
    # coerce to integer
    n <- round(n)
    mode(n) <- "integer"
  } else cli_abort(paste("'n'", error_msg))
  # if p has length 0, make it NULL
  if(!length(p)) p <- NULL
  # when p is a matrix, it must satisfy some conditions
  if(is.matrix(p)) {
    # if matrix has one column, duplicate it
    if(ncol(p) == 1)
      p <- cbind(p, p)
    # stop immediately, if dimensions are violated
    if(ncol(p) != 2)
      cli_abort(paste0("'p' ", error_msg, "; NULL is allowed, too"))
    # check if all values are non-negative and close to integer
    assert_numeric(p, 0, 1, min.len = 2)
  } else
    if(!is.null(p))
      cli_abort(paste0("'p' ", error_msg, "; NULL is allowed, too"))

  len_a <- length(alternative)
  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      tolower(alternative[i]),
      c("two.sided", "less", "greater")
    )
  }

  # method <- match.arg(
  #   tolower(method),
  #   c("pooled", "unpooled", "propdist", "boshloo")
  # )

  qassert(pooled, "B1")

  qassert(exact, "B1")

  if(!exact) qassert(correct, "B1")

  qassert(simple_output, "B1")

  # determine largest number of rows
  len_x <- nrow(x)
  len_n <- nrow(n)
  len_p <- nrow(p)
  len_g <- max(len_x, len_n, len_p, len_a)
  # replicate inputs to same length
  if(len_x < len_g) x <- x[rep_len(seq_len(len_x), len_g), ]
  if(len_n < len_g) n <- n[rep_len(seq_len(len_n), len_g), ]
  if(is.null(p)) p <- matrix(rowSums(x) / rowSums(n), len_g, 2L) else
    if(len_p < len_g) p <- p[rep_len(seq_len(len_p), len_g), ]
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  if(any(x > n))
    cli_abort(paste(
      "Numbers of successes ('x') must not exceed",
      "their respective numbers of trials ('n')."
    ))

  # determine unique parameter sets
  params <- unique(data.frame(n = n, p = p, alternative = alternative))
  n1_u  <- params$n.1
  n2_u  <- params$n.2
  p1_u  <- params$p.1
  p2_u  <- params$p.2
  alt_u <- params$alternative
  len_u <- length(n1_u)

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # begin computations
  for(i in seq_len(len_u)) {
    idx_supp <- which(
      n[, 1] == n1_u[i] & n[, 2] == n2_u[i] &
      p[, 1] == p1_u[i] & p[, 2] == p2_u[i] &
      alternative == alt_u[i]
    )

    # compute test statistics depending on method
    x1 <- rep(0L:n1_u[i], n2_u[i] + 1L) + ifelse(correct && !exact, 0.5, 0)
    x2 <- rep(0L:n2_u[i], rep(n1_u[i] + 1L, n2_u[i] + 1L)) +
      ifelse(correct && !exact, 0.5, 0)
    n1 <- n1_u[i] + ifelse(correct && !exact, 1, 0)
    n2 <- n2_u[i] + ifelse(correct && !exact, 1, 0)
    p1 <- x1 / n1
    p2 <- x2 / n2
    delta <- p1 - p2
    stats <- if(pooled) {
      xx <- x1 + x2
      nn <- n1 + n2
      pp <- xx / nn
      qq <- (nn - xx) / nn
      delta / sqrt(pp * qq * (1/n1 + 1/n2))
    } else {
      q1 <- (n1 - x1) / n1
      q2 <- (n2 - x2) / n2
      delta / sqrt(p1 * q1 / n1 + p2 * q2 / n2)
    }
    stats[delta == 0] <- 0
    stats <- numerical_adjust(stats, FALSE)

    if(exact) {
      # stats <- switch(
      #   EXPR = method,
      #   pooled = {
      #     xx <- x1 + x2
      #     nn <- n1_u[i] + n2_u[i]
      #     pp <- xx / nn
      #     qq <- (nn - xx) / nn
      #     delta / sqrt(pp * qq * (1/n1_u[i] + 1/n2_u[i]))
      #   }
      #   unpooled = {
      #     q1 <- (n1_u[i] - x1) / n1_u[i]
      #     q2 <- (n2_u[i] - x2) / n2_u[i]
      #     delta / sqrt(p1 * q2/n1_u[i] + p2 * q2/n2_u[i])
      #   },
      #   propdist = delta,
      #   boshloo = homogenity_test_pv(
      #     x = cbind(x1, x2),
      #     n = c(n1_u[i], n2_u[i]),
      #     alternative = "two.sided",
      #     ts_method = "minlike",
      #     exact = TRUE,
      #     simple_output = TRUE
      #   )
      # )

      # generate all probabilities under current n and p
      d <- numerical_adjust(
        outer(
          generate_binom_probs(n1_u[i], p1_u[i]),
          generate_binom_probs(n2_u[i], p2_u[i])
        )
      )
      # compute p-value support
      pv_supp <- switch(
        EXPR      = alt_u[i],
        less      = pvalues_by_statistics(stats, d, FALSE),
        greater   = pvalues_by_statistics(stats, d, TRUE),
        two.sided = pvalues_by_statistics(abs(stats), d, TRUE)
      )
    } else {
      if(p1_u[i] * p2_u[i] == 0)
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = c(rep(1, n_u[i] + 1)),
          greater   = c(1, rep(0, n_u[i])),
          two.sided = c(1, rep(0, n_u[i]))
        )
      else if(p1_u[i] * p2_u[i] == 1)
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = c(rep(0, n_u[i]), 1),
          greater   = rep(1, n_u[i] + 1),
          two.sided = c(rep(0, n_u[i]), 1)
        )
      else {
        # compute p-value support
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = pnorm(stats),
          greater   = pnorm(stats, lower.tail = FALSE),
          two.sided = 2 * pnorm(-abs(stats))
        )
      }
    }

    # store results and support
    res[idx_supp] <- pv_supp[x[idx_supp, 2] * (n1_u[i] + 1) + x[idx_supp, 1] + 1]
    if(!simple_output) {
      supports[[i]] <- unique(sort(pv_supp))
      indices[[i]]  <- idx_supp
    }
  }

  return(list(res, supports, indices))
}
