# QR-QURS: mu_tau,t = plogis(X[t,] %*% beta), with fixed tau.
# F(y)=2*sin(pi*y^a/2)/(1+sin(pi*y^a/2)); a=log(C_tau)/log(mu).
# C_p=(2/pi)*asin(p/(2-p)). Sourcing this file defines functions without running the study.
rqq_assert<-function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
rqq_semente<-function(s) {
 rqq_assert(length(s)==1L && is.finite(s) && s>=1 && s<.Machine$integer.max && s==floor(s),
 'Invalid seed.');as.integer(s)
}
rqq_falha<-function(tipo,mensagem,...) stop(structure(c(list(message=mensagem,call=NULL,tipo=tipo),list(...)),
 class=c('erro_rqqurs','error','condition')))
rqq_parametros<-function(mu,tau) {
 rqq_assert(length(tau)==1L && is.finite(tau) && tau>0 && tau<1,'tau must lie in (0, 1).')
 rqq_assert(length(mu)>0 && all(is.finite(mu)) && all(mu>0 & mu<1),'mu must lie in (0, 1).')
 invisible(TRUE)
}
rqq_reciclar<-function(x,n,nome) {
 rqq_assert(length(x)%in%c(1L,n),paste('Invalid length:',nome));rep_len(x,n)
}
rqq_logC<-function(p) log(2/pi)+log(asin(p/(2-p)))
rqq_tag_tau<-function(tau) gsub('.', 'p',sprintf('%.2f',tau),fixed=TRUE)
# Evaluate conditional quantiles at supplied probabilities.
qRQ_QURS<-function(p,mu,tau) {
 rqq_parametros(mu,tau);n<-max(length(p),length(mu));p<-rqq_reciclar(p,n,'p');mu<-rqq_reciclar(mu,n,'mu')
 rqq_assert(all(is.finite(p)) && all(p>=0 & p<=1),'p must lie in [0, 1].')
 ans<-numeric(n);ans[p==1]<-1;i<-which(p>0 & p<1)
 if(length(i)) ans[i]<-exp(rqq_logC(p[i])*log(mu[i])/rqq_logC(tau))
 ans
}
# Evaluate the conditional CDF.
pRQ_QURS<-function(y,mu,tau) {
 rqq_parametros(mu,tau);n<-max(length(y),length(mu));y<-rqq_reciclar(y,n,'y');mu<-rqq_reciclar(mu,n,'mu')
 rqq_assert(all(is.finite(y)),'y must contain finite values.');ans<-numeric(n);ans[y>=1]<-1;i<-which(y>0 & y<1)
 if(length(i)) {a<-rqq_logC(tau)/log(mu[i]);s<-sin((pi/2)*exp(a*log(y[i])));ans[i]<-2*s/(1+s)}
 ans
}
# Evaluate the density or log-density for each conditional quantile.
dRQ_QURS<-function(y,mu,tau,log=FALSE) {
 rqq_parametros(mu,tau);n<-max(length(y),length(mu));y<-rqq_reciclar(y,n,'y');mu<-rqq_reciclar(mu,n,'mu')
 rqq_assert(all(is.finite(y)),'y must contain finite values.');ans<-rep(-Inf,n);i<-which(y>0 & y<1)
 if(length(i)) {
  a<-rqq_logC(tau)/log(mu[i]);ly<-log(y[i]);w<-a*ly;t<-exp(w);z<-(pi/2)*t
  # cos(pi*t/2)=sin(pi*(1-t)/2): stable when t is close to 1.
  cosz<-sin((pi/2)*(-expm1(w)))
  ans[i]<-log(pi)+log(a)+(a-1)*ly+log(cosz)-2*log1p(sin(z))
 }
 if(log) ans else exp(ans)
}
rqq_validar_X<-function(X,n) {
 rqq_assert(is.matrix(X) && is.numeric(X) && nrow(X)==n && ncol(X)>=1 &&
  ncol(X)<n && all(is.finite(X)) && qr(X)$rank==ncol(X),
  'X deve ser matriz numerica finita, de posto completo, com k<n.')
 invisible(TRUE)
}
rqq_nomes<-function(X) paste0('beta_',seq_len(ncol(X))-if(all(X[,1]==1)) 1L else 0L)
# Generate responses conditional on X by inverse transform sampling.
rRQ_QURS<-function(X,beta,tau) {
 rqq_validar_X(X,nrow(X));rqq_assert(length(beta)==ncol(X) && all(is.finite(beta)),'beta must match the design matrix and be finite.')
 mu<-plogis(drop(X%*%beta));rqq_parametros(mu,tau)
 qRQ_QURS(runif(nrow(X)),mu,tau)
}
# Evaluate the negative log-likelihood on the coefficient scale.
rqq_nll<-function(beta,y,X,tau) {
 if(length(beta)!=ncol(X) || any(!is.finite(beta))) return(1e100)
 mu<-plogis(drop(X%*%beta));if(any(mu<=0 | mu>=1)) return(1e100)
 l<-dRQ_QURS(y,mu,tau,log=TRUE);v<--sum(l)
 if(any(!is.finite(l)) || !is.finite(v)) 1e100 else v
}
# Compute the analytic gradient of the negative log-likelihood.
rqq_gradiente<-function(beta,y,X,tau) {
 mu<-plogis(drop(X%*%beta));rqq_parametros(mu,tau)
 L<--log(mu);a<--rqq_logC(tau)/L;w<-a*log(y);z<-(pi/2)*exp(w)
 sn<-sin(z);cs<-sin((pi/2)*(-expm1(w)))
 score_a_vezes_a<-1+w*(1-z*sn/cs-2*z*cs/(1+sn))
 score_eta<-score_a_vezes_a*(1-mu)/L
 g<--drop(crossprod(X,score_eta));rqq_assert(all(is.finite(g)),'Nonfinite gradient.');g
}
# Compute observed information on the coefficient scale.
rqq_hessiana<-function(beta,y,X,tau) {
 p<-length(beta);H<-matrix(NA_real_,p,p)
 for(j in seq_len(p)) {
  h<-.Machine$double.eps^(1/3)*max(1,abs(beta[j]));a<-b<-beta;a[j]<-a[j]+h;b[j]<-b[j]-h
  H[,j]<-(rqq_gradiente(a,y,X,tau)-rqq_gradiente(b,y,X,tau))/(2*h)
 }
 (H+t(H))/2
}
rqq_controle_gosolnp<-function() list(n.restarts=5L,n.sim=200L,sd_inicial=1.5,
 limite_beta=12,expansoes=2L,tol_gradiente=1e-4,tol_rcond=1e-10,margem_preditor=1e-10,
 control=list(trace=0,outer.iter=400L,inner.iter=800L,tol=1e-8,delta=1e-7))
# Fit the model and check the numerical validity of coefficient inference.
estimar_RQ_QURS<-function(y,X,tau,rseed=20260817L,controle=list(),inicial=NULL) {
 if(!requireNamespace('Rsolnp',quietly=TRUE)) stop("Install Rsolnp with install.packages('Rsolnp').")
 rqq_assert(is.numeric(y) && length(y)>0 && all(is.finite(y)) && all(y>0 & y<1),'Data must lie in (0, 1).')
 rqq_validar_X(X,length(y));rqq_parametros(.5,tau);rseed<-rqq_semente(rseed)
 ctrl<-modifyList(rqq_controle_gosolnp(),controle);p<-ncol(X);n<-length(y);nm<-rqq_nomes(X)
 rqq_assert(ctrl$n.restarts>=1 && ctrl$n.sim>=ctrl$n.restarts && ctrl$sd_inicial>0 &&
 ctrl$limite_beta>0 && ctrl$expansoes>=0,'Controles invalidos.')
 # Project a reference constant onto X, keeping the original model matrix.
 if(is.null(inicial)) inicial<-qr.solve(X,rep(qlogis(as.numeric(quantile(y,tau,type=8))),n))
 rqq_assert(length(inicial)==p && all(is.finite(inicial)),'Invalid starting coefficients.')
 inicial<-setNames(as.numeric(inicial),nm);B<-max(ctrl$limite_beta,max(abs(inicial))+2)
 fun<-function(b) rqq_nll(b,y,X,tau)/n;avisos<-character()
 for(expansao in 0:ctrl$expansoes) {
  limite<-B*2^expansao
  seed_local<-as.integer(1+(as.double(rseed)-1+expansao*100000007)%%(.Machine$integer.max-1))
  fit<-withCallingHandlers(Rsolnp::gosolnp(pars=inicial,fun=fun,LB=rep(-limite,p),UB=rep(limite,p),
   distr=rep(2L,p),distr.opt=lapply(inicial,function(m) list(mean=m,sd=ctrl$sd_inicial)),
   n.restarts=ctrl$n.restarts,n.sim=ctrl$n.sim,control=ctrl$control,rseed=seed_local,cluster=NULL),
   warning=function(w) {avisos<<-c(avisos,conditionMessage(w));invokeRestart('muffleWarning')})
  if(length(fit$pars)!=p || any(!is.finite(fit$pars)) || any(abs(fit$pars)>limite))
   rqq_falha('retorno_invalido','gosolnp retornou parametros invalidos.')
  if(!any(abs(fit$pars)>=limite-1e-5)) break
  if(expansao==ctrl$expansoes) rqq_falha('limite_numerico','Coeficiente atingiu limite numerico apos expansoes.')
 }
 if(!isTRUE(fit$convergence==0)) rqq_falha('convergencia',paste('Codigo gosolnp:',fit$convergence),convergencia=fit$convergence)
 beta<-setNames(as.numeric(fit$pars),nm);ll<--rqq_nll(beta,y,X,tau)
 if(!is.finite(ll) || -ll>=1e99) rqq_falha('logverossimilhanca','Nonfinite or penalized fit.')
 mu<-plogis(drop(X%*%beta));eps<-ctrl$margem_preditor
 if(any(mu<=eps | mu>=1-eps)) rqq_falha('limite_numerico','Fitted quantile too close to a numerical boundary.')
 gn<-max(abs(rqq_gradiente(beta,y,X,tau)))/n
 if(gn>ctrl$tol_gradiente) rqq_falha('gradiente','Gradient per observation exceeds the tolerance.',gradiente=gn)
 H<-rqq_hessiana(beta,y,X,tau);dimnames(H)<-list(nm,nm)
 if(any(!is.finite(H))) rqq_falha('hessiana','Nonfinite observed information.')
 ch<-tryCatch(chol(H),error=function(e) NULL)
 if(is.null(ch)) rqq_falha('hessiana','Observed information is not positive definite.')
 rc<-rcond(H);if(!is.finite(rc) || rc<ctrl$tol_rcond) rqq_falha('hessiana','Ill-conditioned observed information.',rcond=rc)
 V<-chol2inv(ch);dimnames(V)<-list(nm,nm);se<-setNames(sqrt(diag(V)),nm)
 if(any(!is.finite(se)) || any(se<=0)) rqq_falha('variancia','Invalid standard errors.')
 lo<-beta-qnorm(.975)*se;hi<-beta+qnorm(.975)*se
 if(any(!is.finite(c(lo,hi))) || any(hi<=lo)) rqq_falha('intervalos','Invalid confidence intervals.')
 list(coef=beta,se=se,lower=lo,upper=hi,vcov=V,hessian=H,loglik=ll,AIC=2*p-2*ll,BIC=log(n)*p-2*ll,
  gradiente_max=gn,rcond=rc,convergencia=fit$convergence,inicial_referencia=inicial,inicial_solver=fit$start.pars,
  rseed=rseed,rseed_efetiva=seed_local,expansoes=expansao,limite_beta=limite,avisos=unique(avisos),
  pacote_Rsolnp=as.character(utils::packageVersion('Rsolnp')))
}
# Predict conditional quantiles with delta-method standard errors and intervals.
rqq_predizer<-function(ajuste,X_novo) {
 rqq_assert(is.matrix(X_novo) && ncol(X_novo)==length(ajuste$coef),'The new design matrix has incompatible dimensions.')
 eta<-drop(X_novo%*%ajuste$coef);v<-rowSums((X_novo%*%ajuste$vcov)*X_novo)
 rqq_assert(all(is.finite(v)) && all(v>0),'Invalid prediction variance.')
 m<-plogis(eta);se_eta<-sqrt(v);se<-m*(1-m)*se_eta;z<-qnorm(.975)
 data.frame(estimativa=m,erro_padrao=se,li_wald=m-z*se,ls_wald=m+z*se,
            li_logit=plogis(eta-z*se_eta),ls_logit=plogis(eta+z*se_eta))
}


