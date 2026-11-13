# QURS raw moments
#
# This script computes E[X^r] for positive orders r under the
# quantile-reparameterized unit rational sine (QURS) distribution.
# No add-on R packages are required.
#
# For a single moment, source this file in a separate environment:
#   qurs <- new.env()
#   sys.source("qurs_moments_corrected.R", envir = qurs)
#   alpha <- qurs$qurs_alpha(mu_tau = 0.50, tau = 0.25)
#   qurs$qurs_moment_integral(r = 1, alpha = alpha)
#
# Run Rscript qurs_moments_corrected.R to reproduce the full comparison. The script
# prints the results and writes qurs_moment_comparison.csv and
# qurs_moment_comparison.tex to the current working directory. The LaTeX
# table uses booktabs commands.
#
# With u = sin((pi/2) * x^alpha) and p = r/alpha, the raw moment is
#   2 * (2/pi)^p * integral_0^1 asin(u)^p / (1 + u)^2 du.
# The Cesaro method replaces (1 + u)^(-2) with a finite weighted sum.
# The comparison uses direct quadrature as its numerical reference.

# Convert the quantile mu_tau and its probability level tau to the positive
# power parameter alpha. Both inputs must be scalars strictly between 0 and 1.
qurs_alpha <- function(mu_tau, tau) {
  stopifnot(length(mu_tau) == 1L, length(tau) == 1L,
            is.finite(mu_tau), is.finite(tau),
            mu_tau > 0, mu_tau < 1, tau > 0, tau < 1)
  A <- (2 / pi) * asin(tau / (2 - tau))
  log(A) / log(mu_tau)
}

# Evaluate the density at x for a given alpha. For finite real inputs, the
# function returns zero at and outside the endpoints of the unit interval.
qurs_pdf <- function(x, alpha) {
  z <- x^alpha
  out <- numeric(length(x))
  ok <- x > 0 & x < 1
  out[ok] <- alpha * pi * x[ok]^(alpha - 1) *
    cos((pi / 2) * z[ok]) / (1 + sin((pi / 2) * z[ok]))^2
  out
}

# Coefficients d_k(p) in [arcsin(u)]^p = u^p sum_k d_k(p) u^(2k).
# a_k are the coefficients of arcsin(u)/u.
qurs_d_coefficients <- function(p, K) {
  stopifnot(length(p) == 1L, p > 0, K >= 0)
  k <- 0:K
  # Use log-gamma arithmetic: direct factorials overflow for K above
  # approximately 170 even though the coefficients themselves are small.
  a <- exp(lgamma(2 * k + 1) - k * log(4) -
             2 * lgamma(k + 1) - log(2 * k + 1))
  d <- numeric(K + 1L)
  d[1L] <- 1
  if (K >= 1L) {
    for (m in 1:K) {
      kk <- 1:m
      d[m + 1L] <- sum((kk * (p + 1) - m) * a[kk + 1L] *
                         d[m - kk + 1L]) / m
    }
  }
  d
}

# Legacy coefficient-based calculation, retained for diagnostics.
# N is the number of outer terms and K is the highest coefficient index.
# For p = r/alpha, the inner sum is
#   H_j = sum_{k=0}^K d_k(p) / (p + j + 2*k + 1).
# Finite K introduces truncation error. Use the kernel implementation to
# evaluate the Cesaro approximation without truncating the inner series.
# Set return_diagnostics = TRUE to inspect the intermediate quantities.
qurs_moment_cesaro_series <- function(r, alpha, N = 4000L, K = 400L,
                                       return_diagnostics = FALSE) {
  stopifnot(r > 0, alpha > 0, N >= 1, K >= 0)
  p <- r / alpha
  d <- qurs_d_coefficients(p, K)
  j <- 0:(N - 1L)
  # Construct the denominator matrix in blocks to limit memory use.
  H <- numeric(N)
  block <- 500L
  for (lo in seq.int(1L, N, by = block)) {
    hi <- min(N, lo + block - 1L)
    jj <- j[lo:hi]
    denom <- outer(jj, 2 * (0:K), "+") + p + 1
    H[lo:hi] <- as.vector((1 / denom) %*% d)
  }
  c <- (j + 1) * (-1)^j * H
  partial <- cumsum(c)
  cesaro <- cumsum(partial) / seq_len(N)
  ans <- 2 * (2 / pi)^p * cesaro[N]
  if (!return_diagnostics) return(ans)
  list(moment = ans, p = p, N = N, K = K, d = d,
       inner_terms = H, outer_terms = c,
       last_partial = partial[N], cesaro = cesaro[N])
}

# Finite Cesaro kernel. If q = -u, then
# K_N(u) = (1/N) sum_{m=0}^{N-1} sum_{j=0}^m (j+1)(-u)^j
#        = sum_{j=0}^{N-1} ((N-j)/N)(j+1)q^j.
# The closed form avoids cancellation from repeatedly accumulating partial
# sums and is stable on 0 <= u <= 1.
qurs_cesaro_kernel <- function(u, N) {
  stopifnot(N >= 1L)
  u <- pmin(pmax(u, 0), 1)
  q <- -u
  n <- 1 - (N + 1) * q^N + N * q^(N + 1)
  den2 <- (1 - q)^2
  den3 <- (1 - q)^3
  np <- ifelse(u == 0, 0, -N * (N + 1) * q^(N - 1) * (1 - q))
  A <- n / den2
  B <- q * (np / den2 + 2 * n / den3)
  A - B / N
}

# Recommended stable Cesaro calculation. The inner quantity H_j(p) is not
# approximated by the d_k series. Instead, the finite Cesaro kernel is
# integrated directly, which isolates the slow outer Cesaro convergence.
qurs_moment_cesaro_kernel <- function(r, alpha, N = 1000L,
                                      subdivisions = 2000L,
                                      rel.tol = 1e-10,
                                      abs.tol = 1e-10) {
  stopifnot(r > 0, alpha > 0, N >= 1L)
  p <- r / alpha
  f <- function(u) {
    2 * (2 / pi)^p * asin(u)^p * qurs_cesaro_kernel(u, N)
  }
  stats::integrate(f, lower = 0, upper = 1,
                   subdivisions = subdivisions,
                   rel.tol = rel.tol, abs.tol = abs.tol)$value
}

# Richardson extrapolation for the leading O(1/N) Cesaro error.
qurs_moment_cesaro_richardson <- function(r, alpha, N = 1000L,
                                          subdivisions = 2000L,
                                          rel.tol = 1e-10,
                                          abs.tol = 1e-10) {
  mN <- qurs_moment_cesaro_kernel(r, alpha, N, subdivisions, rel.tol, abs.tol)
  m2N <- qurs_moment_cesaro_kernel(r, alpha, 2L * N, subdivisions,
                                   rel.tol, abs.tol)
  list(moment = 2 * m2N - mN, moment_N = mN, moment_2N = m2N,
       N = N, N2 = 2L * N)
}

# Default Cesaro interface. N controls the kernel approximation.
# K is retained for compatibility and is unused by the kernel.
# Set return_diagnostics = TRUE to include the inputs and method name.
qurs_moment_cesaro <- function(r, alpha, N = 1000L, K = NULL,
                               return_diagnostics = FALSE) {
  ans <- qurs_moment_cesaro_kernel(r, alpha, N = N)
  if (!return_diagnostics) return(ans)
  list(moment = ans, r = r, alpha = alpha, p = r / alpha,
       N = N, method = "finite Cesaro kernel")
}

# Direct adaptive quadrature of the raw-moment integral on [0, 1].
# The transformation removes the x-scale endpoint factors from the density.
# subdivisions, rel.tol, and abs.tol are passed to stats::integrate().
qurs_moment_integral <- function(r, alpha,
                                 subdivisions = 1000L,
                                 rel.tol = 1e-12, abs.tol = 1e-12) {
  stopifnot(r > 0, alpha > 0)
  p <- r / alpha
  f <- function(u) {
    2 * (2 / pi)^p * (asin(u)^p) / (1 + u)^2
  }
  stats::integrate(f, lower = 0, upper = 1,
                   subdivisions = subdivisions,
                   rel.tol = rel.tol, abs.tol = abs.tol)$value
}

# Double N until successive kernel estimates meet the requested tolerance
# or max_refinements is reached. K is also doubled and recorded for
# compatibility. It does not affect the kernel estimates.
# The result includes the last estimate, convergence flag, and history.
qurs_moment_cesaro_adaptive <- function(r, alpha, N0 = 500L, K0 = 100L,
                                        max_refinements = 6L,
                                        abs.tol = 1e-10, rel.tol = 1e-10) {
  N <- as.integer(N0); K <- as.integer(K0)
  old <- NA_real_
  history <- data.frame(N = integer(), K = integer(), estimate = numeric(),
                        change = numeric())
  for (it in seq_len(max_refinements)) {
    est <- qurs_moment_cesaro(r, alpha, N = N, K = K)
    change <- if (is.na(old)) NA_real_ else abs(est - old)
    history <- rbind(history, data.frame(N = N, K = K,
                                         estimate = est, change = change))
    if (!is.na(change) && change <= abs.tol + rel.tol * abs(est)) {
      return(list(moment = est, N = N, K = K, iterations = it,
                  converged = TRUE, history = history))
    }
    old <- est
    N <- 2L * N
    K <- 2L * K
  }
  list(moment = old, N = N / 2L, K = K / 2L,
       iterations = max_refinements, converged = FALSE, history = history)
}

# Compare raw moments over a data frame with columns tau and mu_tau.
# Each output row represents one parameter pair and one moment order.
# With use.richardson = TRUE, compute the kernel estimates at N and 2N
# and extrapolate. The reported errors use direct quadrature as reference.
# The tolerance arguments apply to the direct quadrature in this comparison.
qurs_compare_moments <- function(cases, moments = 1:3,
                                 N = 4000L, K = 400L,
                                 rel.tol = 1e-12, abs.tol = 1e-12,
                                 use.richardson = TRUE) {
  stopifnot(all(c("tau", "mu_tau") %in% names(cases)))
  out <- list(); z <- 1L
  for (i in seq_len(nrow(cases))) {
    alpha <- qurs_alpha(cases$mu_tau[i], cases$tau[i])
    for (r in moments) {
      series <- qurs_moment_cesaro(r, alpha, N = N, K = K)
      series2 <- if (use.richardson) {
        qurs_moment_cesaro(r, alpha, N = 2L * N, K = K)
      } else NA_real_
      series.richardson <- if (use.richardson) {
        2 * series2 - series
      } else series
      integ <- qurs_moment_integral(r, alpha,
                                    rel.tol = rel.tol, abs.tol = abs.tol)
      out[[z]] <- data.frame(
        tau = cases$tau[i], mu_tau = cases$mu_tau[i], alpha = alpha,
        r = r, moment_cesaro_N = series, moment_cesaro_2N = series2,
        moment_richardson = series.richardson, moment_integral = integ,
        abs_error_N = abs(series - integ),
        abs_error_richardson = abs(series.richardson - integ),
        rel_error_richardson = abs(series.richardson - integ) /
          max(abs(integ), .Machine$double.eps),
        N = N, K = K)
      z <- z + 1L
    }
  }
  do.call(rbind, out)
}

# Example parameter scenarios. Replace or extend this data frame as needed.
cases <- data.frame(
  tau = rep(c(0.10, 0.25, 0.50, 0.75, 0.90), each = 10),
  mu_tau = rep(c(0.05, 0.10, 0.20, 0.30, 0.40,
                 0.50, 0.60, 0.70, 0.80, 0.90),
               times = 5)
)

# Reproduce the comparison for orders 1 through 4 on the example grid
# when run with Rscript or sourced into the global environment.
if (identical(environment(), globalenv())) {
  results <- qurs_compare_moments(cases, moments = 1:4,
                                   N = 1000L, K = 400L,
                                   rel.tol = 1e-12, abs.tol = 1e-12)
  print(results, digits = 10, row.names = FALSE)
  write.csv(results, "qurs_moment_comparison.csv", row.names = FALSE)

  # Compact LaTeX table for the Supplementary Material.
  latex <- capture.output({
    cat("% Generated by qurs_moments_corrected.R\n")
    cat("\\begin{tabular}{rrrrrrrrr}\n")
    cat("\\toprule\n")
    cat("$\\tau$ & $\\mu_\\tau$ & $r$ & $\\alpha$ & $M_N$ & $M_{2N}$ & $M_R$ & Integration & $|\\Delta_R|$\\\\\n")
    cat("\\midrule\n")
    for (i in seq_len(nrow(results))) {
      x <- results[i, ]
      cat(sprintf("%.2f & %.2f & %d & %.6f & %.10f & %.10f & %.10f & %.10f & %.3e\\\\\n",
                  x$tau, x$mu_tau, x$r, x$alpha, x$moment_cesaro_N,
                  x$moment_cesaro_2N, x$moment_richardson,
                  x$moment_integral, x$abs_error_richardson))
    }
    cat("\\bottomrule\n\\end{tabular}\n")
  })
  writeLines(latex, "qurs_moment_comparison.tex")
}

