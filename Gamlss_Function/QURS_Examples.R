# Reproducible QURS examples
# From the directory containing both files, run:
#   Rscript QURS_examples_fixed.R
# Install dependencies once if needed:
#   install.packages(c("gamlss", "gamlss.dist"))

library(gamlss)
source("QURS_GAMLSS.R")

# Example 1: Fit a QURS distribution with a common 25th percentile.
set.seed(2026)
tau_distribution <- 0.25
true_quantile <- 0.60
distribution_data <- data.frame(
  y = rQURS(400, mu = true_quantile, tau = tau_distribution)
)

distribution_fit <- gamlss(
  y ~ 1, family = QURS(tau = tau_distribution),
  data = distribution_data, control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(distribution_fit$converged))
cat("Example 1: Common 25th percentile\n")
print(c(
  true_quantile = true_quantile,
  estimated_quantile = unname(fitted(distribution_fit, what = "mu")[1])
))

# Example 2: Fit a QR-QURS model for the conditional 75th percentile.
# The logit of the conditional quantile is linear in x and group.
set.seed(2027)
tau_regression <- 0.75
n <- 600
regression_data <- data.frame(
  x = stats::runif(n, -1, 1),
  group = stats::rbinom(n, size = 1, prob = 0.5)
)
true_beta <- c("(Intercept)" = -0.3, x = 0.8, group = 0.4)
conditional_quantile <- with(
  regression_data, stats::plogis(true_beta[1] + true_beta[2] * x +
                                  true_beta[3] * group)
)
regression_data$y <- rQURS(n, mu = conditional_quantile, tau = tau_regression)

regression_fit <- gamlss(
  y ~ x + group, family = QURS(tau = tau_regression),
  data = regression_data, control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(regression_fit$converged))
cat("\nExample 2: Coefficients on the logit scale\n")
print(cbind(true = true_beta, estimated = coef(regression_fit, what = "mu")))

# Fitted conditional quantiles are on the original response scale.
cat("\nFirst six fitted conditional 75th percentiles\n")
print(head(fitted(regression_fit, what = "mu")))
