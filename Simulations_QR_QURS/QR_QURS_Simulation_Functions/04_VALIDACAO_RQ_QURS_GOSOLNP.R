# Validate the distribution, derivatives, optimizer, and checkpoint recovery on small runs.
.rqq_pasta<-local({
 a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
 fs<-Filter(function(x)!is.null(x),lapply(sys.frames(),function(x)x$ofile))
 f<-if(length(fs)) tail(fs,1)[[1]] else if(length(a)) sub('^--file=','',a[1]) else NULL
 if(is.null(f)) getwd() else dirname(normalizePath(f,mustWork=TRUE))
})
source(file.path(.rqq_pasta,'01_FUNCOES_RQ_QURS_GOSOLNP.R'),encoding='UTF-8')
source(file.path(.rqq_pasta,'02_MOTOR_MONTE_CARLO_RQ_QURS.R'),encoding='UTF-8')
if(!requireNamespace('Rsolnp',quietly=TRUE)) stop("Install Rsolnp with install.packages('Rsolnp').")
out<-file.path(.rqq_pasta,'VALIDACAO_RQ_QURS_GOSOLNP');dir.create(out,showWarnings=FALSE)
checks<-character();ok<-function(msg){checks<<-c(checks,msg);message('OK: ',msg)}
taus<-c(.1,.25,.5,.75,.9);err<-0
for(mu in c(.05,.1,.5,.9,.95)) for(tau in taus) {
 p<-c(.001,.01,.1,.25,.5,.75,.9,.99,.999)
 q<-qRQ_QURS(p,mu,tau);err<-max(err,max(abs(pRQ_QURS(q,mu,tau)-p)))
 stopifnot(abs(qRQ_QURS(tau,mu,tau)-mu)<1e-12)
 v<-integrate(function(y)dRQ_QURS(y,mu,tau),0,1,subdivisions=2000L,rel.tol=1e-7)
 stopifnot(abs(v$value-1)<1e-5)
 y<-c(.1,.3,.5,.7,.9);a<-rqq_logC(tau)/log(mu)
 closed<-pi*a*y^(a-1)*cos(pi*y^a/2)/(1+sin(pi*y^a/2))^2
 stopifnot(max(abs(dRQ_QURS(y,mu,tau)-closed))<1e-8)
}
stopifnot(err<1e-9);ok(sprintf('25 parameter settings: quantile, CDF, normalization, and closed-form density; inversion error %.3g',err))
cs<-rqq_cenarios();Xbase<-rqq_desenhos_fixos(c(50L,100L,200L,400L));g<-rqq_grade(cs)
stopifnot(nrow(g)==80L,nrow(rqq_grade(rqq_cenarios(FALSE)))==60L)
stopifnot(nrow(rqq_selecionar_grade(g,cenario_inicial='R3'))==40L,
 nrow(rqq_selecionar_grade(g,cenarios_selecionados=c('R1','R3')))==40L)
for(c in cs) {
 a<-rqq_matrizes_cenario(c,Xbase[['50']]);b<-rqq_matrizes_cenario(c,Xbase[['400']])
 stopifnot(identical(a,b[1:50,,drop=FALSE]))
}
stopifnot(identical(rqq_nomes(rqq_matrizes_cenario(cs$R1,Xbase[['50']])),paste0('beta_',0:2)),
 identical(rqq_nomes(rqq_matrizes_cenario(cs$R2,Xbase[['50']])),paste0('beta_',1:4)),
 identical(rqq_nomes(rqq_matrizes_cenario(cs$R3,Xbase[['50']])),paste0('beta_',1:6)))
stopifnot(identical(vapply(taus,rqq_tag_tau,''),c('0p10','0p25','0p50','0p75','0p90')))
ok('Verified scenarios, intercepts, fixed designs, selection, and distinct labels for tau.')
ctrl<-list(n.restarts=3L,n.sim=100L);cases<-list();records<-list();maxg<-0
for(i in seq_along(cs)) for(tau in taus) {
 c<-cs[[i]];X<-rqq_matrizes_cenario(c,Xbase[['100']]);beta<-c$beta
 set.seed(18000L+100L*i+round(100*tau));y<-rRQ_QURS(X,beta,tau)
 h<-1e-5;ng<-vapply(seq_along(beta),function(j){
  a<-b<-beta;a[j]<-a[j]+h;b[j]<-b[j]-h
  (rqq_nll(a,y,X,tau)-rqq_nll(b,y,X,tau))/(2*h)
 },0.0)
 ge<-max(abs(ng-rqq_gradiente(beta,y,X,tau)));maxg<-max(maxg,ge);stopifnot(ge<1e-4)
 f<-tryCatch(estimar_RQ_QURS(y,X,tau,rseed=40000L+100L*i+round(100*tau),controle=ctrl),error=function(e)e)
 records[[length(records)+1L]]<-data.frame(cenario=c$cenario,tau=tau,n=100,
  status=if(inherits(f,'error'))'REJEITADO' else 'OK',erro=if(inherits(f,'error'))conditionMessage(f) else '',
  erro_gradiente=ge,gradiente_ajuste=if(inherits(f,'error'))NA_real_ else f$gradiente_max)
 cases[[length(cases)+1L]]<-list(fit=f,y=y,X=X,tau=tau)
}
rqq_csv(rqq_linhas(records),file.path(out,'Ajustes_20_Casos.csv'))
ok(sprintf('Gradient checks in 20 cases: maximum error %.3g',maxg))
valid<-Filter(function(a)!inherits(a$fit,'error'),cases);stopifnot(length(valid)>0)
a<-valid[[1]];f<-a$fit;f2<-estimar_RQ_QURS(a$y,a$X,a$tau,rseed=f$rseed,controle=ctrl)
stopifnot(identical(f$coef,f2$coef),identical(f$se,f2$se))
ok('The same optimizer seed reproduces coefficients and standard errors.')
if(requireNamespace('numDeriv',quietly=TRUE)) {
 H<-numDeriv::hessian(function(b)rqq_nll(b,a$y,a$X,a$tau),f$coef)
 e<-max(abs(H-f$hessian))/max(1,max(abs(H)));stopifnot(e<1e-4)
 ok(sprintf('Observed Hessian compared with numDeriv: relative error %.3g',e))
}
cfg<-g[g$cenario=='R2' & g$tau==.5 & g$n==100,,drop=FALSE]
gc<-0L;ec<-0L
genfail<-function(...){gc<<-gc+1L;if(gc==1L)stop('FALHA_GERACAO_TESTE');rRQ_QURS(...)}
estfail<-function(...){ec<<-ec+1L;if(ec==1L)rqq_falha('hessiana','FALHA_HESSIANA_TESTE');estimar_RQ_QURS(...)}
tmp<-tempfile('rqqurs_validacao_')
st<-rqq_executar_configuracao(cfg,cs$R2,Xbase[['100']],tmp,B=3L,controle=ctrl,checkpoint_cada=1L,
 max_novas_tentativas=2L,gerador=genfail,estimador=estfail)
stopifnot(length(st$validos)==0L,length(st$tentativas)==2L)
st<-rqq_executar_configuracao(cfg,cs$R2,Xbase[['100']],tmp,B=3L,controle=ctrl,checkpoint_cada=1L)
stopifnot(length(st$validos)==3L)
at<-rqq_linhas(st$tentativas);stopifnot(any(at$erro=='FALHA_GERACAO_TESTE'),
 any(at$erro=='FALHA_HESSIANA_TESTE' & at$tipo=='hessiana'))
st2<-rqq_executar_configuracao(cfg,cs$R2,Xbase[['100']],tmp,B=3L,controle=ctrl)
stopifnot(identical(st$validos,st2$validos),identical(st$tentativas,st2$tentativas))
s<-read.csv(file.path(tmp,'Resumo_Monte_Carlo_RQ_QURS.csv'));co<-read.csv(file.path(tmp,'Controle_Replicas.csv'))
stopifnot(nrow(s)==4L,co$completo,co$replicas_validas==3L,co$tentativas-co$replicas_descartadas==3L)
for(j in seq_along(cs$R2$beta)) {
 rr<-s[s$parametro==paste0('beta_',j),,drop=FALSE]
 v<-vapply(st$validos,function(f)f$coef[j],0.0);lo<-vapply(st$validos,function(f)f$lower[j],0.0);hi<-vapply(st$validos,function(f)f$upper[j],0.0)
 stopifnot(abs(rr$media-mean(v))<1e-12,abs(rr$eqm-mean((v-cs$R2$beta[j])^2))<1e-12,
 rr$cobertura_wald==mean(lo<=cs$R2$beta[j] & cs$R2$beta[j]<=hi))
}
bad<-try(rqq_executar_configuracao(cfg,cs$R2,Xbase[['100']],tmp,B=4L,controle=ctrl),silent=TRUE)
stopifnot(inherits(bad,'try-error'));unlink(tmp,recursive=TRUE)
ok('Verified failure logging, replacement, summary metrics, checkpoint recovery, and rejection of incompatible checkpoints.')
capture.output(sessionInfo(),file=file.path(out,'SessionInfo.txt'))
writeLines(c('QR-QURS VALIDATION',checks,paste('Valid fits across the 20 cases:',length(valid)),
 'This validation does not run the full study.'),file.path(out,'Relatorio_Validacao.txt'))
message('Validation complete: ',out)


