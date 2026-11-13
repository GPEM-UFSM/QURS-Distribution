# Quantile Unit Rational Sine (QURS) distribution for GAMLSS
#
# mu is the quantile at the fixed probability level tau. Both are in (0, 1).
# The power parameter is alpha = log((2/pi)*asin(tau/(2-tau))) / log(mu).
# See QURS_examples_fixed.R for an intercept-only fit and a QR-QURS regression.

qurs_log_density <- expression(log(pi)+log((log((2/pi)*asin(tau/(2-tau)))/log(mu)))+
                    ((log((2/pi)*asin(tau/(2-tau)))/log(mu))-1)*log(y)+
                    log(cos((pi*y^(log((2/pi)*asin(tau/(2-tau)))/log(mu)))/2))-
                    2*log(1+sin((pi*y^(log((2/pi)*asin(tau/(2-tau)))/log(mu)))/2)))

# Symbolic derivatives with respect to mu on the response scale.
qurs_score <- D(qurs_log_density, "mu")
qurs_curvature <- D(qurs_score, "mu")

qurs_validate_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau <= 0 || tau >= 1) {
    stop("tau must be a single number strictly between 0 and 1.")
  }
}

qurs_family_alpha <- function(mu, tau) {
  qurs_validate_tau(tau)
  if (!is.numeric(mu) || !length(mu) ||
      any(!is.finite(mu) | mu <= 0 | mu >= 1)) {
    stop("mu must contain finite quantiles strictly between 0 and 1.")
  }
  log((2 / pi) * asin(tau / (2 - tau))) / log(mu)
}

# Construct a family with one modeled parameter, mu. The quantile level
# is captured by the family and stays fixed throughout fitting.
# The default logit link maps (0, 1) onto the real line.
QURS <- function(mu.link = "logit", tau = 0.50) {
  qurs_validate_tau(tau)
  force(tau)
  link <- gamlss.dist::checklink(
    "mu.link", "QURS", substitute(mu.link),
    c("logit", "probit", "cloglog", "identity", "log", "own")
  )

  # GAMLSS evaluates the deviance body in its own fitting environment.
  # Embed each evaluator as a function object in the callback body so
  # the fixed quantile level remains available after that extraction.
  bind_callback <- function(evaluator) {
    callback <- function(y, mu, ...) NULL
    body(callback) <- as.call(list(evaluator, quote(y), quote(mu)))
    callback
  }

  score <- bind_callback(function(y, mu) {
    eval(qurs_score, list(y = y, mu = mu, tau = tau))
  })
  curvature <- bind_callback(function(y, mu) {
    # Retain the original negative curvature floor for fitting.
    value <- eval(qurs_curvature, list(y = y, mu = mu, tau = tau))
    ifelse(value < -1e-15, value, -1e-15)
  })
  deviance <- bind_callback(function(y, mu) {
    -2 * dQURS(y, mu, tau, log = TRUE)
  })
  quantile_residuals <- bind_callback(function(y, mu) {
    stats::qnorm(pQURS(q = y, mu = mu, tau = tau))
  })

  structure(
    list(
      family = c("QURS", "Quantile Unit Rational Sine"),
      parameters = list(mu = TRUE),
      nopar = 1,
      type = "Continuous",
      mu.link = as.character(substitute(mu.link)),
      mu.linkfun = link$linkfun,
      mu.linkinv = link$linkinv,
      mu.dr = link$mu.eta,
      dldm = score,
      d2ldm2 = curvature,
      G.dev.incr = deviance,
      rqres = as.expression(body(quantile_residuals)),
      mu.initial = expression(mu <- (y + mean(y)) / 2),
      mu.valid = function(mu) all(is.finite(mu) & mu > 0 & mu < 1),
      y.valid = function(y) all(is.finite(y) & y > 0 & y < 1)
    ),
    class = c("gamlss.family", "family")
  )
}

# Density on the open unit interval. Set log = TRUE for the log-density.
dQURS <- function(y, mu, tau = 0.50, log = FALSE) {
  qurs_family_alpha(mu, tau)
  if (any(!is.finite(y) | y <= 0 | y >= 1)) {
    stop("y must contain finite observations strictly between 0 and 1.")
  }
  log_density <- eval(qurs_log_density)
  if (log) log_density else exp(log_density)
}

# Cumulative distribution function, including the support endpoints.
# Use lower.tail = FALSE for survival probabilities and log.p = TRUE
# for probabilities on the logarithmic scale.
pQURS <- function(q, mu, tau = 0.50, lower.tail = TRUE, log.p = FALSE) {
  alpha <- qurs_family_alpha(mu, tau)
  theta <- (pi / 2) * pmin(pmax(q, 0), 1)^alpha
  sine <- sin(theta)
  probability <- if (lower.tail) {
    2 * sine / (1 + sine)
  } else {
    # Equivalent to 1 - F(q), with less cancellation near one.
    ifelse(theta == pi / 2, 0, (cos(theta) / (1 + sine))^2)
  }
  if (log.p) log(probability) else probability
}

# Quantile function for probabilities in [0, 1]. In particular,
# qQURS(tau, mu, tau) equals mu up to floating-point rounding.
qQURS <- function(p, mu, tau = 0.50) {
  alpha <- qurs_family_alpha(mu, tau)
  if (any(!is.finite(p) | p < 0 | p > 1)) {
    stop("p must contain finite probabilities between 0 and 1.")
  }
  ((2 / pi) * asin(p / (2 - p)))^(1 / alpha)
}

# Generate n observations by inverse transform sampling.
# Supply either a common mu or one conditional quantile per observation.
rQURS <- function(n, mu, tau = 0.50) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n)) {
    stop("n must be a nonnegative integer.")
  }
  if (!length(mu) %in% c(1L, n) || !length(mu)) {
    stop("mu must be a scalar or contain one quantile per observation.")
  }
  qQURS(stats::runif(n), mu, tau)
}
