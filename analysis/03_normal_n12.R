# The small-sample setting of the original study: N(theta, 1), n = 12, theta = 1,
# one-sided alpha = 0.05. Shows what changes when every test has the same size.
# Run from the repository root: Rscript analysis/03_normal_n12.R
source("R/signed_rank.R")
source("R/power.R")
dir.create("results", showWarnings = FALSE)

set.seed(1984)
n <- 12
theta <- 1
alpha <- 0.05
reps <- 200000
samples <- matrix(rnorm(n * reps, mean = theta), nrow = n)

# Non-randomised versions: reject when T >= smallest k with P(T >= k) <= alpha
conservative <- function(pmf) {
  upper <- rev(cumsum(rev(pmf)))
  k <- which(upper <= alpha)[1] - 1
  list(k = k, size = upper[k + 1])
}
sr <- conservative(signed_rank_pmf(n))
sg <- conservative(dbinom(0:n, n, 0.5))
p_pos <- pnorm(theta)

out <- data.frame(
  test = c("z (known variance)", "t", "signed rank", "signed rank, randomised",
           "sign", "sign, randomised"),
  size = c(alpha, alpha, sr$size, alpha, sg$size, alpha),
  power = c(
    pnorm(qnorm(1 - alpha) - sqrt(n) * theta, lower.tail = FALSE),
    pt(qt(1 - alpha, n - 1), n - 1, ncp = sqrt(n) * theta, lower.tail = FALSE),
    mean(signed_rank_statistic(samples) >= sr$k),
    mean(reject_signed_rank(samples, alpha)),
    pbinom(sg$k - 1, n, p_pos, lower.tail = FALSE),
    mean(reject_sign(samples, alpha))
  ),
  method = c("exact", "exact", "simulation", "simulation", "exact", "simulation")
)
write.csv(out, "results/normal_n12.csv", row.names = FALSE)
print(out, digits = 4, row.names = FALSE)
cat("\nCritical values: signed rank T+ >=", sr$k, "| sign S >=", sg$k, "\n")
