# One-sided tests of H0: theta = 0 against H1: theta > 0, applied column-wise to
# a matrix of samples (n rows, one column per replication). Each function returns
# the rejection probability per column: 0, 1, or the randomisation weight on the
# boundary of a discrete test, so that every discrete test has size exactly alpha.

# Critical value k and weight gamma such that P(T > k) + gamma * P(T = k) = alpha.
randomised_cutoff <- function(pmf, alpha) {
  tail_above <- rev(cumsum(rev(pmf))) - pmf      # P(T > t), t = 0, 1, ...
  i <- which(tail_above <= alpha)[1]
  list(k = i - 1, gamma = (alpha - tail_above[i]) / pmf[i])
}

critical_function <- function(stat, cut) {
  (stat > cut$k) + cut$gamma * (stat == cut$k)
}

signed_rank_statistic <- function(samples) {
  ranks <- apply(abs(samples), 2, rank)
  colSums(ranks * (samples > 0))
}

reject_signed_rank <- function(samples, alpha = 0.05) {
  cut <- randomised_cutoff(signed_rank_pmf(nrow(samples)), alpha)
  critical_function(signed_rank_statistic(samples), cut)
}

reject_sign <- function(samples, alpha = 0.05) {
  n <- nrow(samples)
  cut <- randomised_cutoff(dbinom(0:n, n, 0.5), alpha)
  critical_function(colSums(samples > 0), cut)
}

reject_t <- function(samples, alpha = 0.05) {
  n <- nrow(samples)
  m <- colMeans(samples)
  s <- sqrt(colSums(sweep(samples, 2, m)^2) / (n - 1))
  as.numeric(sqrt(n) * m / s > qt(1 - alpha, n - 1))
}

# Symmetric error distributions centred at 0. All have unit variance except the
# Cauchy, which is scaled to share the interquartile range of the standard normal.
cn_sd <- sqrt(0.9 * 1 + 0.1 * 25)

error_distributions <- list(
  normal = list(
    r = function(m) rnorm(m),
    d = function(x) dnorm(x)),
  uniform = list(
    r = function(m) runif(m, -sqrt(3), sqrt(3)),
    d = function(x) dunif(x, -sqrt(3), sqrt(3))),
  laplace = list(
    r = function(m) (rexp(m) - rexp(m)) / sqrt(2),
    d = function(x) exp(-sqrt(2) * abs(x)) / sqrt(2)),
  student_t3 = list(
    r = function(m) rt(m, 3) / sqrt(3),
    d = function(x) sqrt(3) * dt(sqrt(3) * x, 3)),
  contaminated_normal = list(
    r = function(m) rnorm(m, sd = ifelse(runif(m) < 0.1, 5, 1)) / cn_sd,
    d = function(x) cn_sd * (0.9 * dnorm(cn_sd * x) + 0.1 * dnorm(cn_sd * x, sd = 5))),
  cauchy = list(
    r = function(m) rcauchy(m, scale = qnorm(0.75)),
    d = function(x) dcauchy(x, scale = qnorm(0.75)))
)

# Pitman asymptotic relative efficiencies with respect to the t test, for a
# symmetric density f with variance sigma2:
#   signed rank: 12 sigma2 (integral of f^2)^2      sign: 4 sigma2 f(0)^2
efficiency_vs_t <- function(d, sigma2 = 1) {
  f2 <- function(x) d(x)^2
  limit <- if (d(10) == 0) sqrt(3) else Inf      # compact support: uniform
  integral <- 2 * integrate(f2, 0, limit, rel.tol = 1e-10)$value
  c(signed_rank = 12 * sigma2 * integral^2, sign = 4 * sigma2 * d(0)^2)
}
