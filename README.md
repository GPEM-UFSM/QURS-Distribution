# QURS Distribution and QR-QURS Regression Model

This repository provides R code and Monte Carlo replication files for the
Quantile Unit Rational Sine (QURS) distribution and its quantile regression
model (QR-QURS). The distribution is defined for continuous responses on
the open unit interval $(0,1)$. Its parameter $\mu_\tau$ is the quantile
at a fixed probability level $\tau$. In QR-QURS regression, covariates
enter through a link function for the conditional quantile.

## Repository contents

| Directory | Contents |
|:--|:--|
| [Gamlss_Function](Gamlss_Function/) | QURS family for GAMLSS and examples of fitting an intercept-only distribution and a conditional quantile regression. |
| [Moments](Moments/) | R code for computing raw QURS moments and comparing adaptive quadrature, finite Cesàro approximations, and Richardson extrapolation. |
| [Simulations_QURS](Simulations_QURS/) | R scripts for the QURS simulation and post-processing, with 100 archived replication files for scenarios C1–C5. |
| [Simulations_QR_QURS](Simulations_QR_QURS/) | R scripts for QR-QURS simulation, numerical validation, and post-processing, with 60 archived replication files for scenarios R1–R3. |
| [Shiny_App](Shiny_App/) | Interactive plots of the density, cumulative distribution, survival, and hazard functions, along with Bowley skewness and Moors kurtosis. |
| [Launch Shiny App](https://gpem-ufsm.github.io/QURS-Distribution/) | Run the interactive QURS application directly in your browser. |

## Getting started

The scripts and application use R. Open the README in the relevant
directory for package requirements and examples.

- **Fit a model:** From [Gamlss_Function](Gamlss_Function/), run
  `Rscript QURS_Examples.R`.
- **Compute moments:** From [Moments](Moments/), run
  `Rscript QURS_Moments.R` to generate the CSV and LaTeX comparison tables.
- **Explore the distribution:** From the repository root, run
  `shiny::runApp("Shiny_App")` in R after installing the app dependencies.

## Monte Carlo replication files

The [QURS archive](Simulations_QURS/QURS_Simulations/) and the
[QR-QURS archive](Simulations_QR_QURS/QR-QURS_Simulations/) contain one
`.rds` file per configuration. Each configuration combines a scenario,
a sample size, and a quantile level and contains 5,000 accepted
replications. From the repository root, use `readRDS()` to inspect an
archived configuration:

```r
qurs <- readRDS("Simulations_QURS/QURS_Simulations/C1_tau_01_n_050.rds")
qr_qurs <- readRDS("Simulations_QR_QURS/QR-QURS_Simulations/R1_tau_01_n_050.rds")
```

Both examples select the first scenario at $\tau=0.10$ and $n=50$.
The simulation READMEs describe the scenario parameters and how to rerun
the studies and generate tables and figures.

