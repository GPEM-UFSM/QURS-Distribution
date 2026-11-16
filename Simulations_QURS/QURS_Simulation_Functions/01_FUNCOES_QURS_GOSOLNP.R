# QURS distribution and maximum likelihood estimation using Rsolnp::gosolnp.
# The distribution is parameterized by its quantile mu_tau at level tau.
# alpha = log{(2/pi) asin[tau/(2-tau)]}/log(mu_tau) > 0.
# F(y) = 2 sin(pi*y^alpha/2) / {1 + sin(pi*y^alpha/2)}.
# Pass tau explicitly to each function. Sourcing this file only defines functions.

.validar_prob_QURS <- function(x, nome) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop(nome, ' must be a scalar in (0, 1).', call. = FALSE)
}
.log_base_QURS <- function(p) log((2/pi)*asin(p/(2-p)))
# Convert the quantile parameter to the power parameter alpha.
alpha_QURS <- function(mu_tau, tau) {
  .validar_prob_QURS(mu_tau, 'mu_tau'); .validar_prob_QURS(tau, 'tau')
  .log_base_QURS(tau)/log(mu_tau)
}
.log_dQURS_eta <- function(y, eta) {
  alpha <- exp(eta)
  if (!is.finite(alpha) || alpha <= 0) return(rep(-Inf, length(y)))
  ly <- log(y)
  t <- -alpha*ly
  z <- (pi/2)*exp(-t)
  # cos(z) = sin(pi/2-z); expm1 avoids cancellation when y^alpha is close to 1.
  delta <- (pi/2)*(-expm1(-t))
  log(pi) + eta + (alpha-1)*ly + log(sin(delta)) - 2*log1p(sin(z))
}
# Evaluate the density or log-density at a vector of points.
dQURS <- function(x, mu_tau, tau, log = FALSE) {
  alpha <- alpha_QURS(mu_tau, tau)
  if (!is.numeric(x)) stop('x must be numeric.')
  ans <- rep(-Inf, length(x)); ans[is.na(x)] <- NA_real_
  j <- which(!is.na(x) & x>0 & x<1)
  ans[j] <- .log_dQURS_eta(x[j], log(alpha))
  # One-sided density limits at 0 and 1.
  j0 <- which(!is.na(x) & x==0)
  ans[j0] <- if (alpha<1) Inf else if (alpha==1) log(pi) else -Inf
  if (log) ans else exp(ans)
}
# Evaluate the CDF or survival function, optionally on the log scale.
pQURS <- function(q, mu_tau, tau, lower.tail = TRUE, log.p = FALSE) {
  alpha <- alpha_QURS(mu_tau, tau)
  if (!is.numeric(q)) stop('q must be numeric.')
  ans <- rep(NA_real_, length(q))
  ans[!is.na(q) & q<=0] <- if (lower.tail) -Inf else 0
  ans[!is.na(q) & q>=1] <- if (lower.tail) 0 else -Inf
  j <- which(!is.na(q) & q>0 & q<1)
  t <- -alpha*log(q[j]); z <- (pi/2)*exp(-t)
  delta <- (pi/2)*(-expm1(-t))
  if (lower.tail) {
    # sin(z)/z -> 1 as z approaches zero. Compute log(z) without underflow.
    sinc <- rep(1, length(z)); k <- which(z>1e-8)
    sinc[k] <- sin(z[k])/z[k]
    ans[j] <- pmin(0, log(2)+log(pi/2)-t+log(sinc)-log1p(sin(z)))
  } else {
    # 1-F = (1-sin(z))/(1+sin(z)) = [cos(z)/(1+sin(z))]^2.
    ans[j] <- 2*(log(sin(delta))-log1p(sin(z)))
  }
  if (log.p) ans else exp(ans)
}
# Evaluate the quantile function for probabilities in [0, 1].
qQURS <- function(p, mu_tau, tau) {
  alpha <- alpha_QURS(mu_tau, tau)
  if (!is.numeric(p) || any(p<0 | p>1, na.rm=TRUE)) stop('p must lie in [0, 1].')
  exp(.log_base_QURS(p)/alpha)
}
# Generate a sample with a common mu_tau by inverse transform sampling.
rQURS <- function(n, mu_tau, tau) {
  if (!is.numeric(n) || length(n)!=1L || !is.finite(n) || n<0 || n!=floor(n))
    stop('n must be a nonnegative integer.')
  qQURS(runif(n), mu_tau, tau)
}
# Sum the log-density over the sample.
logLik_QURS <- function(dados, mu_tau, tau) {
  if (!is.numeric(dados) || !length(dados) || any(!is.finite(dados)) ||
      any(dados<=0 | dados>=1)) stop('Data must be finite and strictly between 0 and 1.')
  sum(dQURS(dados, mu_tau, tau, log=TRUE))
}
# Compute the sample score with respect to eta = log(alpha).
.score_eta_QURS <- function(eta, dados) {
  t <- -exp(eta)*log(dados)
  z <- (pi/2)*exp(-t); delta <- (pi/2)*(-expm1(-t))
  cz <- sin(delta); sz <- sin(z)
  sum(1-t + (t*z)*(sz/cz + 2*cz/(1+sz)))
}
# Derive a separate optimizer seed from the sample seed.
# The estimator optimizes mu_tau directly inside the unit interval.
.semente_otimizacao_QURS <- function(semente_amostra) as.integer(1+
  (as.double(semente_amostra)+1000000000) %% (.Machine$integer.max-1))

# Fit the model and return the estimate, observed information, and diagnostics.
# dados is a numeric sample on (0, 1); tau is the fixed quantile level.
estimar_QURS <- function(dados,tau,controle=list(),rseed=20261016L,
    n.restarts=5L,n.sim=100L,margem_numerica=1e-8) {
  if(!requireNamespace('Rsolnp',quietly=TRUE))
    stop("Install Rsolnp with install.packages('Rsolnp') before running this script.")
  .validar_prob_QURS(tau,'tau')
  if(!is.numeric(dados) || !length(dados) || any(!is.finite(dados)) ||
     any(dados<=0 | dados>=1)) stop('Data must be finite and strictly between 0 and 1.')
  for(x in list(n.restarts,n.sim))
    if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<1 || x!=floor(x))
      stop('n.restarts and n.sim must be positive integers.')
  if(n.sim<n.restarts) stop('n.sim must be at least n.restarts.')
  if(!is.numeric(rseed) || length(rseed)!=1L || !is.finite(rseed) ||
     rseed<0 || rseed!=floor(rseed) || rseed>=.Machine$integer.max)
    stop('rseed must be a valid integer seed.')
  if(length(margem_numerica)!=1L || !is.finite(margem_numerica) ||
     margem_numerica<=0 || margem_numerica>=.1) stop('Invalid numerical boundary margin.')
  controle <- modifyList(list(trace=0,outer.iter=400L,inner.iter=800L,
    tol=1e-8,delta=1e-7),controle)
  lower <- margem_numerica; upper <- 1-margem_numerica
  nll_mu <- function(theta) {
    if(length(theta)!=1L || !is.finite(theta) || theta<=0 || theta>=1) return(1e100)
    eta <- log(alpha_QURS(unname(theta),tau))
    v <- -mean(.log_dQURS_eta(dados,eta))
    if(!is.finite(v)) 1e100 else v
  }
  # Let gosolnp generate its own random starting points.
  sol <- Rsolnp::gosolnp(pars=NULL,fun=nll_mu,LB=lower,UB=upper,distr=1L,
    n.restarts=n.restarts,n.sim=n.sim,rseed=as.integer(rseed),cluster=NULL,
    control=controle)
  if(!is.list(sol) || length(sol$pars)!=1L || length(sol$convergence)!=1L)
    stop('Invalid result from Rsolnp::gosolnp.')
  mu <- unname(as.numeric(sol$pars)); convergencia <- as.integer(sol$convergence)
  if(!is.finite(mu) || mu<=0 || mu>=1 || mu<lower || mu>upper)
    stop('gosolnp returned a parameter outside the numerical bounds.')
  alpha <- alpha_QURS(mu,tau); eta <- log(alpha)
  ll <- sum(.log_dQURS_eta(dados,eta))
  # Recompute observed information for log(alpha) and apply the delta
  # method to obtain standard errors for mu_tau and alpha.
  nll_eta <- function(e) -sum(.log_dQURS_eta(dados,e))
  h <- unname(stats::optimHess(eta,nll_eta)[1,1])
  hp <- is.finite(h) && h>0
  var_eta <- if(hp) 1/h else NA_real_
  dmu <- -mu*log(mu)
  se_mu <- dmu*sqrt(var_eta); se_alpha <- alpha*sqrt(var_eta)
  score <- .score_eta_QURS(eta,dados)
  mensagens <- character()
  if(convergencia!=0L) mensagens <- c(mensagens,paste('gosolnp: codigo',convergencia))
  if(!is.finite(ll) || nll_mu(mu)>=1e99) {
    convergencia <- 90L;mensagens <- c(mensagens,'Nonfinite likelihood or penalized region.')
  }
  if(!is.finite(score) || abs(score)>1e-4*max(1,length(dados))) {
    convergencia <- 91L;mensagens <- c(mensagens,'Score exceeds the verification tolerance.')
  }
  no_limite <- mu<=lower+1e-10 || mu>=upper-1e-10
  regular <- convergencia==0L && !no_limite && hp && is.finite(se_mu) && se_mu>0
  if(no_limite) mensagens <- c(mensagens,'LIMITE_NUMERICO_MU: regular inference is unavailable.')
  structure(list(estimativas=c(mu_tau=mu),erro_padrao=c(mu_tau=se_mu),
    alpha_estimado=alpha,se_alpha=se_alpha,log_alpha_estimado=eta,
    vcov=matrix(se_mu^2,1,1,dimnames=list('mu_tau','mu_tau')),
    var_log_alpha=var_eta,informacao_log_alpha=h,escore_log_alpha=score,
    hessiana_positiva=hp,inferencia_regular=regular,convergencia=convergencia,
    logLik=ll,tau=tau,n=length(dados),mensagem=paste(mensagens,collapse=' | '),
    metodo='Rsolnp::gosolnp',semente_otimizacao=as.integer(rseed),
    n_restarts=n.restarts,n_sim=n.sim,limites=c(lower,upper),
    versao_Rsolnp=as.character(utils::packageVersion('Rsolnp')),
    hessiana_gosolnp=sol$hessian,otimizacao=sol),class='ajuste_QURS')
}

# Return a two-column confidence interval for mu_tau. The logit method
# constructs the interval on the logit scale and transforms its endpoints.
confint_QURS <- function(ajuste, nivel=0.95, metodo=c('wald','logit')) {
  metodo <- match.arg(metodo)
  .validar_prob_QURS(nivel,'nivel')
  if (!isTRUE(ajuste$inferencia_regular)) stop('Regular inference is unavailable.')
  mu <- unname(ajuste$estimativas['mu_tau']); se <- unname(ajuste$erro_padrao['mu_tau'])
  z <- qnorm((1+nivel)/2)
  ic <- if (metodo=='wald') mu+c(-1,1)*z*se else
    plogis(qlogis(mu)+c(-1,1)*z*se/(mu*(1-mu)))
  matrix(ic,1,2,dimnames=list('mu_tau',c('inferior','superior')))
}


