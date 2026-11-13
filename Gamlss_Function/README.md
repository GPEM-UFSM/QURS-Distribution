# QURS and QR-QURS with GAMLSS

`QURS_GAMLSS.R` defines the Quantile Unit Rational Sine (QURS) family
for the R package `gamlss`. Use an intercept-only model to fit a common
quantile, or include covariates to fit a QR-QURS regression.

## Setup

Install the required packages once:

```r
install.packages(c("gamlss", "gamlss.dist"))
```

Run the examples from the directory containing both R scripts:

```sh
Rscript QURS_Examples.R
```

To work interactively from that directory:

```r
library(gamlss)
source("QURS_GAMLSS.R")
```

## Family and parameterization

The response must be finite and strictly between zero and one.
In `QURS(mu.link = "logit", tau = 0.50)`, `tau` is a fixed scalar
in (0, 1). The modeled parameter `mu` is the quantile at probability
level `tau`. For instance, setting `tau = 0.50` models the conditional
median. Under the default logit link, `logit(mu) = X %*% beta`.
The fitted values for `mu` are on the original response scale,
and the regression coefficients are on the logit scale.

## Example 1. Fit a QURS distribution

Generate 400 observations with a true 25th percentile of 0.60 and
estimate that percentile with an intercept-only model.

```r
set.seed(2026)
tau_distribution <- 0.25
distribution_data <- data.frame(
  y = rQURS(400, mu = 0.60, tau = tau_distribution)
)

distribution_fit <- gamlss(
  y ~ 1,
  family = QURS(tau = tau_distribution),
  data = distribution_data,
  control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(distribution_fit$converged))
fitted(distribution_fit, what = "mu")[1]
```

The fitted value is the estimated 25th percentile and is common to
all observations.

## Example 2. Fit a QR-QURS regression

Generate 600 observations whose conditional 75th percentile depends on
a continuous covariate `x` and a binary indicator `group`.

```r
set.seed(2027)
tau_regression <- 0.75
n <- 600
regression_data <- data.frame(
  x = runif(n, -1, 1),
  group = rbinom(n, size = 1, prob = 0.5)
)

mu <- with(regression_data, plogis(-0.3 + 0.8 * x + 0.4 * group))
regression_data$y <- rQURS(n, mu = mu, tau = tau_regression)

regression_fit <- gamlss(
  y ~ x + group,
  family = QURS(tau = tau_regression),
  data = regression_data,
  control = gamlss.control(trace = FALSE)
)
stopifnot(isTRUE(regression_fit$converged))
coef(regression_fit, what = "mu")
head(fitted(regression_fit, what = "mu"))
```

Holding the other covariate fixed, a one-unit increase in a predictor
changes the logit of the conditional 75th percentile by its
coefficient. The fitted values are estimated conditional 75th
percentiles. Use `residuals(regression_fit)` to obtain the normalized
quantile residuals for diagnostic plots.

## Distribution functions

| Function | Purpose |
|:--|:--|
| `dQURS(y, mu, tau, log = FALSE)` | Density for observations in (0, 1), or log-density when `log = TRUE`. |
| `pQURS(q, mu, tau, lower.tail = TRUE, log.p = FALSE)` | Cumulative probability or, with `lower.tail = FALSE`, upper-tail probability. |
| `qQURS(p, mu, tau)` | Quantile at probability `p` in [0, 1]. In particular, `qQURS(tau, mu, tau)` returns `mu` up to rounding. |
| `rQURS(n, mu, tau)` | Generate `n` observations with one common `mu` or one conditional `mu` per observation. |

These helpers default to `tau = 0.50`. Pass the same `tau` to the data
generator and fitted family when reproducing a simulation.
