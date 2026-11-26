# Run from the repository root: Rscript tests/run_tests.R
source("R/signed_rank.R")
source("R/power.R")

near <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))

# 1. Null distribution against R's own implementation
for (n in c(1, 2, 5, 10, 25, 50)) {
  t <- 0:(n * (n + 1) / 2)
  stopifnot(near(dsr(t, n), dsignrank(t, n)), near(psr(t, n), psignrank(t, n)))
}

# 2. Structural properties
for (n in 1:15) {
  p <- signed_rank_pmf(n)
  t <- 0:(n * (n + 1) / 2)
  stopifnot(
    near(sum(p), 1),
    near(p, rev(p)),
    near(sum(t * p), n * (n + 1) / 4),
    near(sum((t - n * (n + 1) / 4)^2 * p), n * (n + 1) * (2 * n + 1) / 24)
  )
}

# 3. Out-of-support arguments
stopifnot(dsr(-1, 4) == 0, dsr(11, 4) == 0, dsr(2.5, 4) == 0,
          psr(-1, 4) == 0, psr(99, 4) == 1)

# 4. p-values against stats::wilcox.test, all alternatives
set.seed(2024)
for (i in 1:500) {
  n <- sample(2:30, 1)
  x <- rnorm(n, mean = runif(1, -1, 1))
  mu <- runif(1, -1, 1)
  for (alt in c("two.sided", "less", "greater")) {
    ours <- signed_rank_test(x, mu = mu, alternative = alt)
    base <- wilcox.test(x, mu = mu, alternative = alt, exact = TRUE)
    stopifnot(near(unname(ours$statistic), unname(base$statistic)),
              near(ours$p.value, base$p.value))
  }
}

# 5. Extreme samples: T+ = 0 and T+ = n(n+1)/2
low <- -(1:6)
high <- 1:6
stopifnot(
  near(signed_rank_test(low, alternative = "greater")$p.value, 1),
  near(signed_rank_test(low, alternative = "less")$p.value, 2^-6),
  near(signed_rank_test(high, alternative = "greater")$p.value, 2^-6),
  near(signed_rank_test(high, alternative = "less")$p.value, 1)
)

# 6. Ties are refused, zeros are dropped
stopifnot(inherits(try(signed_rank_test(c(1, -1, 2)), silent = TRUE), "try-error"))
stopifnot(near(
  suppressMessages(signed_rank_test(c(0, 1, 2, 3)))$p.value,
  signed_rank_test(c(1, 2, 3))$p.value
))

# 7. Randomised cut-offs give size exactly alpha
for (n in c(5, 12, 20, 41)) {
  for (alpha in c(0.01, 0.05, 0.1)) {
    for (pmf in list(signed_rank_pmf(n), dbinom(0:n, n, 0.5))) {
      cut <- randomised_cutoff(pmf, alpha)
      support <- seq_along(pmf) - 1
      stopifnot(cut$gamma >= 0, cut$gamma <= 1,
                near(sum(critical_function(support, cut) * pmf), alpha))
    }
  }
}

# 8. Efficiencies against their closed forms
eff <- function(name) efficiency_vs_t(error_distributions[[name]]$d)
stopifnot(
  near(eff("normal"), c(signed_rank = 3 / pi, sign = 2 / pi), 1e-6),
  near(eff("laplace"), c(signed_rank = 1.5, sign = 2), 1e-6),
  near(eff("uniform"), c(signed_rank = 1, sign = 1 / 3), 1e-6)
)

cat("All tests passed.\n")
