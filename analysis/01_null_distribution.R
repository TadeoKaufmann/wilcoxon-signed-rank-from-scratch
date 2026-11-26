# Exact null distribution of T+ against its normal approximation.
# Run from the repository root: Rscript analysis/01_null_distribution.R
source("R/signed_rank.R")
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

moments <- function(n) c(mean = n * (n + 1) / 4,
                         sd = sqrt(n * (n + 1) * (2 * n + 1) / 24))

png("figures/null_distribution.png", width = 1500, height = 480, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
for (n in c(4, 10, 20)) {
  t <- 0:(n * (n + 1) / 2)
  m <- moments(n)
  p <- signed_rank_pmf(n)
  plot(t, p, type = "h", lwd = if (n < 20) 4 else 2,
       ylim = c(0, 1.05 * max(p, dnorm(0, sd = m["sd"]))),
       col = "steelblue", xlab = expression(T^"+"), ylab = "Probability",
       main = paste("n =", n))
  curve(dnorm(x, m["mean"], m["sd"]), add = TRUE, col = "firebrick", lwd = 2, n = 500)
}
invisible(dev.off())

# Largest CDF error of the normal approximation, with and without continuity correction
ns <- c(4, 6, 8, 10, 15, 20, 30, 50, 100)
errors <- do.call(rbind, lapply(ns, function(n) {
  t <- 0:(n * (n + 1) / 2)
  m <- moments(n)
  exact <- psr(t, n)
  data.frame(
    n = n,
    plain = max(abs(exact - pnorm(t, m["mean"], m["sd"]))),
    continuity_corrected = max(abs(exact - pnorm(t + 0.5, m["mean"], m["sd"])))
  )
}))
write.csv(errors, "results/normal_approximation_error.csv", row.names = FALSE)
print(errors, digits = 3, row.names = FALSE)

png("figures/normal_approximation_error.png", width = 900, height = 560, res = 130)
par(mar = c(4, 4.5, 1, 1))
matplot(errors$n, errors[, -1], type = "b", log = "xy", pch = c(16, 17), lty = 1,
        col = c("steelblue", "firebrick"), xlab = "n (log scale)",
        ylab = "Max |exact CDF - normal CDF| (log scale)")
legend("topright", c("Plain", "Continuity corrected"), pch = c(16, 17), lty = 1,
       col = c("steelblue", "firebrick"), bty = "n")
invisible(dev.off())
