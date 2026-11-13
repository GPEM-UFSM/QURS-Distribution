# QURS Shiny App

This R Shiny app explores the Quantile Unit Rational Sine (QURS)
distribution on the open unit interval. Use it to visualize
quantile-based shape measures and compare the density, distribution,
survival, and hazard curves for selected parameters.

## Run locally

Install the required R packages once:

```r
install.packages(c("shiny", "bslib", "ggplot2", "purrr", "viridis", "dplyr"))
```

From the directory containing [app.R](app.R), start the app:

```r
shiny::runApp(".")
```

If your working directory is the repository root, use
`shiny::runApp("Shiny_App")`.

## Parameters and controls

The parameter `mu_tau` is the quantile at probability level `tau`, so
`F(mu_tau) = tau`. For example, `tau = 0.50` makes `mu_tau` the median.

| Tab | Controls | Output |
|:--|:--|:--|
| **Skewness and Kurtosis** | Upper limits for `mu_tau` and `tau`. Each grid starts at 0.01. | Bowley skewness and Moors kurtosis as functions of `mu_tau`, with curves grouped by `tau`. |
| **Distribution** | A fixed `mu_tau`, a fixed `tau), and the upper endpoint of the plotted response range, `y`. | Probability density, cumulative distribution, survival, and hazard for that parameter pair. |

### Example

In **Distribution**, set `tau = 0.50`, `mu_tau = 0.70`, and
`y = 0.99`. The cumulative distribution curve passes through
`F(0.70) = 0.50`, and the survival probability at the same response
value is also 0.50. Change `mu_tau` with `tau` fixed to see how all
four curves change.

In **Skewness and Kurtosis**, set both upper limits to 0.90 to
compare the shape measures across a wider range of parameter values.
Bowley skewness is based on quartiles, and Moors kurtosis is based on
octiles.
