# Wilcoxon signed-rank test, implemented from first principles (base R only).

# Null probabilities P(T+ = t), t = 0, ..., n(n+1)/2.
# Under H0, T+ = sum_{j=1}^n j * W_j with W_j iid Bernoulli(1/2), so the pmf is
# built by adding one rank at a time: O(n^3) instead of the O(2^n) recursion.
signed_rank_pmf <- function(n) {
  stopifnot(length(n) == 1, n >= 0, n == round(n))
  p <- 1
  for (k in seq_len(n)) {
    p <- (c(p, numeric(k)) + c(numeric(k), p)) / 2
  }
  p
}

dsr <- function(t, n) {
  p <- signed_rank_pmf(n)
  out <- numeric(length(t))
  ok <- t == round(t) & t >= 0 & t <= length(p) - 1
  out[ok] <- p[t[ok] + 1]
  out
}

psr <- function(t, n) {
  cdf <- cumsum(signed_rank_pmf(n))
  top <- length(cdf) - 1
  t <- floor(t)
  out <- numeric(length(t))
  inside <- t >= 0 & t <= top
  out[inside] <- cdf[t[inside] + 1]
  out[t > top] <- 1
  pmin(out, 1)
}

# Exact test for the centre of a symmetric continuous distribution.
signed_rank_test <- function(x, mu = 0,
                             alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  data_name <- deparse(substitute(x))
  stopifnot(is.numeric(x), length(mu) == 1, is.finite(mu))

  y <- x[is.finite(x)] - mu
  if (any(y == 0)) {
    message("Dropping ", sum(y == 0), " observation(s) equal to mu.")
    y <- y[y != 0]
  }
  n <- length(y)
  if (n == 0) stop("No observations left to test.")
  if (anyDuplicated(abs(y))) {
    stop("Tied absolute values: the exact null distribution does not apply. ",
         "See stats::wilcox.test for a normal approximation with tie correction.")
  }

  t_plus <- sum(rank(abs(y))[y > 0])
  lower <- psr(t_plus, n)                        # P(T+ <= t)
  upper <- sum(dsr(t_plus:(n * (n + 1) / 2), n)) # P(T+ >= t)

  p_value <- switch(alternative,
    less = lower,
    greater = min(1, upper),
    two.sided = min(1, 2 * min(lower, upper))
  )

  structure(list(
    statistic = c(V = t_plus),
    p.value = p_value,
    null.value = c(location = mu),
    alternative = alternative,
    method = "Exact Wilcoxon signed-rank test",
    data.name = data_name
  ), class = "htest")
}
