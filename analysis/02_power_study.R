# Size and power of the t, signed-rank and sign tests across error distributions.
# Run from the repository root: Rscript analysis/02_power_study.R
source("R/signed_rank.R")
source("R/power.R")
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

set.seed(1984)
n <- 20
reps <- 40000
alpha <- 0.05
thetas <- seq(0, 1, by = 0.1)
tests <- list(t = reject_t, signed_rank = reject_signed_rank, sign = reject_sign)

power <- do.call(rbind, lapply(names(error_distributions), function(dist) {
  # Common random numbers: the same errors are shifted by every theta
  errors <- matrix(error_distributions[[dist]]$r(n * reps), nrow = n)
  do.call(rbind, lapply(thetas, function(theta) {
    samples <- errors + theta
    do.call(rbind, lapply(names(tests), function(test) {
      phi <- tests[[test]](samples, alpha)
      data.frame(distribution = dist, theta = theta, test = test,
                 power = mean(phi), mc_se = sd(phi) / sqrt(reps))
    }))
  }))
}))
write.csv(power, "results/power.csv", row.names = FALSE)

wide <- function(theta) {
  s <- power[abs(power$theta - theta) < 1e-9, ]
  tab <- with(s, tapply(power, list(distribution, test), identity))
  tab[names(error_distributions), names(tests)]
}
cat("\nEmpirical size (theta = 0, nominal 0.05)\n"); print(round(wide(0), 4))
cat("\nPower at theta = 0.5\n"); print(round(wide(0.5), 4))

efficiency <- t(sapply(setdiff(names(error_distributions), "cauchy"), function(dist) {
  efficiency_vs_t(error_distributions[[dist]]$d)
}))
write.csv(data.frame(distribution = rownames(efficiency), efficiency),
          "results/efficiency_vs_t.csv", row.names = FALSE)
cat("\nAsymptotic relative efficiency vs t test\n"); print(round(efficiency, 3))

colours <- c(t = "grey30", signed_rank = "steelblue", sign = "darkorange")
png("figures/power_curves.png", width = 1500, height = 950, res = 130)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
for (dist in names(error_distributions)) {
  plot(NA, xlim = range(thetas), ylim = c(0, 1), xlab = expression(theta),
       ylab = "Power", main = gsub("_", " ", dist))
  abline(h = alpha, lty = 3, col = "grey60")
  for (test in names(tests)) {
    s <- power[power$distribution == dist & power$test == test, ]
    lines(s$theta, s$power, col = colours[test], lwd = 2)
  }
  if (dist == "normal") {
    legend("bottomright", c("t", "signed rank", "sign"), col = colours,
           lwd = 2, bty = "n")
  }
}
invisible(dev.off())
