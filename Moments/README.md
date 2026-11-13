# Raw moments of the QURS distribution

[QURS_Moments.R](QURS_Moments.R) computes raw moments $m_r = \mathrm{E}[X^r]$
of positive order $r$ for the Quantile Unit Rational Sine (QURS)
distribution on $(0,1)$. It uses standard R packages only.

## Parameters

- `tau` is a fixed probability level in $(0,1)$.
- `mu_tau` is the corresponding quantile, also in $(0,1)$.
- `alpha` is the positive power parameter obtained with
  `qurs_alpha(mu_tau, tau)`. The individual moment functions take
  `alpha` as input.

With $p=r/\alpha$ and $u=\sin(\pi x^\alpha/2)$, the moment integral is

$$
m_r = 2(2/\pi)^p \int_0^1 \frac{[\arcsin(u)]^p}{(1+u)^2}\,du.
$$

## Example 1. Mean and variance

From the directory containing the script, load the functions into a
separate environment to avoid triggering the full comparison and file exports.

```r
qurs <- new.env()
sys.source("QURS_Moments.R", envir = qurs)

alpha <- qurs$qurs_alpha(mu_tau = 0.50, tau = 0.25)
m1 <- qurs$qurs_moment_integral(r = 1, alpha = alpha)
m2 <- qurs$qurs_moment_integral(r = 2, alpha = alpha)
c(mean = m1, variance = m2 - m1^2)
```

The first raw moment is the mean. The variance is the second raw moment
minus the square of the first.

## Example 2. Compare numerical methods

```r
cases <- data.frame(tau = c(0.25, 0.75), mu_tau = c(0.50, 0.70))
comparison <- qurs$qurs_compare_moments(cases, moments = 1:4, N = 1000L)

comparison[c("tau", "mu_tau", "r", "moment_integral",
             "moment_cesaro_N", "moment_richardson",
             "abs_error_richardson")]
```

Each row represents a parameter pair and a moment order. The error
columns measure differences relative to direct quadrature.

| Function | Method |
|:--|:--|
| `qurs_moment_integral(r, alpha)` | Adaptive quadrature of the transformed integral. This is the reference used in the comparisons. |
| `qurs_moment_cesaro(r, alpha, N)` | Quadrature with the finite Cesàro kernel. `N` controls the number of terms. Use `return_diagnostics = TRUE` for a list with the estimate and settings. |
| `qurs_moment_cesaro_richardson(r, alpha, N)` | Computes the kernel estimates at `N` and `2N` and returns both along with the extrapolation `2 * M_2N - M_N`. |
| `qurs_moment_cesaro_series(r, alpha, N, K)` | Evaluates the coefficient-based series. `K` is the largest coefficient index retained and introduces a separate truncation error. |
| `qurs_moment_cesaro_adaptive(r, alpha)` | Doubles `N` until successive kernel estimates meet the tolerance or the refinement limit is reached. Check `converged` and `history` in the result. |

The kernel methods accept or report `K` in some interfaces for
compatibility. Changing it affects the coefficient-based series only.
The quadrature tolerances control integration error. The choice of
`N` determines the finite-kernel approximation, and agreement between
successive estimates is a stopping criterion rather than an error bound.

## Run the full comparison

```sh
Rscript QURS_Moments.R
```

This evaluates moments of orders 1–4 for 50 pairs formed by crossing
`tau = 0.10, 0.25, 0.50, 0.75, 0.90` with
`mu_tau = 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90`.
The script prints the 200 rows and writes the following files to the
current working directory:

- `qurs_moment_comparison.csv` contains the estimates and discrepancies.
- `qurs_moment_comparison.tex` contains a LaTeX table fragment. Load
  `booktabs` in the parent document.

Running the script again replaces these outputs.
