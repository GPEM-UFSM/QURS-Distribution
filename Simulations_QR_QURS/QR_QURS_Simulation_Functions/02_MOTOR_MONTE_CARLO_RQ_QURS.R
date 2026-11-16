# Monte Carlo engine. Sourcing this file does not run the study.
# Source 01_FUNCOES_RQ_QURS_GOSOLNP.R first.

# Define coefficient vectors and covariate structures.
rqq_cenarios <- function(incluir_adicional=TRUE) {
 out<-list(
  R1=list(cenario='R1',beta=c(.5,.3,.2),intercepto=TRUE,rho=0,
          tipo='com_intercepto',descricao='R1 solicitado: beta_0,beta_1,beta_2'),
  R2=list(cenario='R2',beta=c(.05,-.27,-.81,.05),intercepto=FALSE,rho=0,
          tipo='sem_intercepto',descricao='R2 solicitado: beta_1 a beta_4, sem intercepto'),
  R3=list(cenario='R3',beta=c(.045,.03,-.27,-.05,-.81,.05),intercepto=FALSE,rho=0,
          tipo='sem_intercepto',descricao='R3 solicitado: beta_1 a beta_6, sem intercepto'))
 if(incluir_adicional) out$R4<-list(cenario='R4',beta=c(.5,.3,.2),intercepto=TRUE,rho=.8,
          tipo='covariaveis_correlacionadas',descricao='Comparacao com R1; copula Gaussiana com rho latente .8')
 out
}

# Cross the scenarios with sample sizes and quantile levels.
rqq_grade <- function(cenarios,taus=c(.1,.25,.5,.75,.9),ns=c(50L,100L,200L,400L)) {
  rqq_assert(length(cenarios)>0 && !anyDuplicated(names(cenarios)),'Cenarios invalidos.')
  rqq_assert(length(taus)>0 && !anyDuplicated(taus) && all(taus>0 & taus<1),'taus invalidos.')
  rqq_assert(!anyDuplicated(vapply(taus,rqq_tag_tau,'')),'Niveis tau colidem nos nomes de pasta.')
  rqq_assert(length(ns)>0 && !anyDuplicated(ns) && all(ns==floor(ns)) && all(ns>6),'ns invalidos.')
  g<-expand.grid(n=ns,tau=taus,cenario=names(cenarios),stringsAsFactors=FALSE)
  g$id_configuracao<-seq_len(nrow(g));g[,c('id_configuracao','cenario','tau','n')]
}
rqq_selecionar_grade <- function(grade,cenario_inicial=NULL,cenarios_selecionados=NULL) {
  disponiveis<-unique(grade$cenario)
  rqq_assert(is.null(cenario_inicial)||is.null(cenarios_selecionados),
    'Use cenario_inicial OU cenarios_selecionados.')
  if(!is.null(cenario_inicial)) {
    rqq_assert(length(cenario_inicial)==1L && cenario_inicial%in%disponiveis,'Cenario inicial invalido.')
    cenarios_selecionados<-disponiveis[seq.int(match(cenario_inicial,disponiveis),length(disponiveis))]
  }
  if(is.null(cenarios_selecionados)) return(grade)
  rqq_assert(length(cenarios_selecionados)>0 && all(cenarios_selecionados%in%disponiveis),
    'Cenario selecionado inexistente.')
  grade[grade$cenario%in%cenarios_selecionados,,drop=FALSE]
}
# Draw one base design and take nested row subsets.
rqq_desenhos_fixos <- function(ns,semente=220260817L) {
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(rqq_semente(semente))
 nmax<-max(ns);base<-matrix(runif(nmax*6L),nmax,6L)
 colnames(base)<-paste0('x',1:6)
 out<-lapply(ns,function(n) base[seq_len(n),,drop=FALSE]);names(out)<-as.character(ns);out
}
# Construct the design matrix for the requested scenario.
rqq_matrizes_cenario <- function(cenario,base) {
 p<-length(cenario$beta)-as.integer(cenario$intercepto)
 rqq_assert(p>=1 && p<=ncol(base),'Quantidade de covariaveis incompativel.')
 X<-base[,seq_len(p),drop=FALSE]
 if(cenario$rho!=0) {
  rqq_assert(p>=2 && abs(cenario$rho)<1,'Correlacao invalida.')
  X[,2]<-pnorm(cenario$rho*qnorm(X[,1])+sqrt(1-cenario$rho^2)*qnorm(X[,2]))
 }
 X<-if(cenario$intercepto) cbind(intercepto=1,X) else X
 rqq_validar_X(X,nrow(X));X
}
rqq_verdade <- function(cenario) cenario$beta

rqq_csv <- function(x,arquivo) utils::write.csv(x,arquivo,row.names=FALSE,na='',fileEncoding='UTF-8')
rqq_salvar_rds <- function(x,arquivo) {
  tmp<-paste0(arquivo,'.tmp');saveRDS(x,tmp,version=2)
  if(file.exists(arquivo)) rqq_assert(file.copy(arquivo,paste0(arquivo,'.bak'),overwrite=TRUE),
    'Falha ao preservar checkpoint anterior.')
  if(!file.rename(tmp,arquivo)) {
    rqq_assert(file.copy(tmp,arquivo,overwrite=TRUE),'Falha ao salvar checkpoint.');unlink(tmp)
  }
  invisible(arquivo)
}
rqq_ler_checkpoint <- function(arquivo) {
  a<-tryCatch(readRDS(arquivo),error=function(e) NULL)
  if(is.null(a) && file.exists(paste0(arquivo,'.bak'))) {
    a<-tryCatch(readRDS(paste0(arquivo,'.bak')),error=function(e) NULL)
    if(!is.null(a)) warning('Checkpoint principal ilegivel; retomando copia anterior.')
  }
  rqq_assert(!is.null(a),'Checkpoint ilegivel: revise os arquivos antes de continuar.');a
}
rqq_linhas <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
rqq_erro_df <- function() data.frame(tentativa=integer(),semente_amostra=integer(),
  semente_otimizacao=integer(),valida=logical(),tipo=character(),erro=character(),
  avisos=character(),tempo=double(),convergencia=integer(),gradiente_max=double(),
  rcond=double(),stringsAsFactors=FALSE)
rqq_estatisticas <- function(est,se,lo,hi,verdade) {
  verdade<-unname(verdade)
  m<-length(est);rqq_assert(m>=1L,'Sem estimativas para resumir.')
  centro<-est-mean(est);m2<-mean(centro^2)
  c(verdadeiro=verdade,replicas=m,media=mean(est),vies=mean(est)-verdade,
    vies_relativo_percentual=if(verdade!=0) 100*(mean(est)-verdade)/verdade else NA_real_,
    desvio_padrao=if(m>1) sd(est) else NA_real_,erro_padrao_medio=mean(se),
    eqm=mean((est-verdade)^2),raiz_eqm=sqrt(mean((est-verdade)^2)),
    assimetria=if(m2>0) mean(centro^3)/m2^1.5 else NA_real_,
    curtose=if(m2>0) mean(centro^4)/m2^2 else NA_real_,
    cobertura_wald=mean(lo<=verdade & verdade<=hi),
    nao_cobertura_wald=mean(!(lo<=verdade & verdade<=hi)),amplitude_media_wald=mean(hi-lo))
}
rqq_resumir_coeficientes <- function(validos,verdade,config,tipo) {
  if(!length(validos)) return(data.frame())
  out<-lapply(seq_along(verdade),function(j) {
    est<-vapply(validos,function(x) x$coef[j],0.0)
    se<-vapply(validos,function(x) x$se[j],0.0)
    lo<-vapply(validos,function(x) x$lower[j],0.0);hi<-vapply(validos,function(x) x$upper[j],0.0)
    cbind(config,tipo_modelo=tipo,parametro=names(verdade)[j],
      escala='coeficiente_regressao',
      tipo_intervalo='wald_coeficiente',
      as.data.frame(as.list(rqq_estatisticas(est,se,lo,hi,verdade[j]))),row.names=NULL)
  });do.call(rbind,out)
}
rqq_pontos_referencia <- function(cenario) {
 p<-length(cenario$beta)-as.integer(cenario$intercepto)
 X<-matrix(rep(c(.25,.5,.75),p),3L,p);colnames(X)<-paste0('x',seq_len(p))
 list(X=if(cenario$intercepto) cbind(intercepto=1,X) else X,
      nomes=c('x025','x050','x075'))
}
rqq_resumir_predicoes <- function(validos,cenario,config) {
 if(!length(validos)) return(data.frame())
 ref<-rqq_pontos_referencia(cenario);pred<-lapply(validos,function(f) rqq_predizer(f,ref$X))
 true<-plogis(drop(ref$X%*%cenario$beta));out<-list()
 for(j in seq_len(3L)) {
  e<-vapply(pred,function(x) x$estimativa[j],0.0);se<-vapply(pred,function(x) x$erro_padrao[j],0.0)
  lw<-vapply(pred,function(x) x$li_wald[j],0.0);uw<-vapply(pred,function(x) x$ls_wald[j],0.0)
  ll<-vapply(pred,function(x) x$li_logit[j],0.0);ul<-vapply(pred,function(x) x$ls_logit[j],0.0)
  row<-cbind(config,tipo_modelo=cenario$tipo,parametro='mu_tau',ponto=ref$nomes[j],
    covariaveis=paste(ref$X[j,],collapse=';'),as.data.frame(as.list(rqq_estatisticas(e,se,lw,uw,true[j]))))
  row$cobertura_logit<-mean(ll<=true[j] & true[j]<=ul);row$amplitude_media_logit<-mean(ul-ll)
  out[[j]]<-row
 };do.call(rbind,out)
}

rqq_validos_csv <- function(validos) {
  if(!length(validos)) return(data.frame())
  rqq_linhas(lapply(seq_along(validos),function(i) {
    x<-validos[[i]]
    row<-data.frame(replica_valida=i,tentativa=x$tentativa,semente_amostra=x$semente_amostra,
      semente_otimizacao=x$rseed,semente_solver_efetiva=x$rseed_efetiva,
      tempo=x$tempo,loglik=x$loglik,AIC=x$AIC,BIC=x$BIC,gradiente_max=x$gradiente_max,
      rcond=x$rcond,expansoes=x$expansoes,limite_beta=x$limite_beta,
      avisos=paste(x$avisos,collapse=' | '),stringsAsFactors=FALSE)
    for(nm in names(x$coef)) {
      row[[paste0('estimativa_',nm)]]<-x$coef[nm];row[[paste0('ep_',nm)]]<-x$se[nm]
      row[[paste0('li_wald_',nm)]]<-x$lower[nm];row[[paste0('ls_wald_',nm)]]<-x$upper[nm]
      row[[paste0('inicial_referencia_',nm)]]<-x$inicial_referencia[nm]
      st<-as.numeric(x$inicial_solver)
      row[[paste0('inicial_solver_',nm)]]<-if(length(st)==length(x$coef)) st[match(nm,names(x$coef))] else NA_real_
    };row
  }))
}
rqq_exportar_configuracao <- function(estado,pasta,config,cenario,X,verdade) {
  dir.create(pasta,recursive=TRUE,showWarnings=FALSE)
  val<-estado$validos;at<-if(length(estado$tentativas)) rqq_linhas(estado$tentativas) else rqq_erro_df()
  descartadas<-at[!at$valida,,drop=FALSE];N<-nrow(at);B<-estado$assinatura$B
  completo<-length(val)==B
  controle<-cbind(config,tipo_modelo=cenario$tipo,replicas_solicitadas=B,
    replicas_validas=length(val),tentativas=N,replicas_descartadas=nrow(descartadas),
    taxa_descarte=if(N) nrow(descartadas)/N else 0,completo=completo,
    tempo_total_tentativas=sum(at$tempo),tempo_medio_replica_valida=if(length(val))
      mean(vapply(val,function(v) v$tempo,0.0)) else NA_real_)
  rqq_csv(controle,file.path(pasta,'Controle_Replicas.csv'))
  rqq_csv(at,file.path(pasta,'Registro_Tentativas.csv'))
  rqq_csv(descartadas,file.path(pasta,'Replicas_Descartadas.csv'))
  rqq_csv(rqq_validos_csv(val),file.path(pasta,'Replicas_Validas.csv'))
  saveRDS(list(configuracao=config,cenario=cenario,X=X,
    verdade=verdade,assinatura=estado$assinatura,replicas=val,completo=completo),
    file.path(pasta,'Replicas_Validas.rds'),version=2)
  rqq_csv(rqq_resumir_coeficientes(val,verdade,config,cenario$tipo),
    file.path(pasta,'Resumo_Monte_Carlo_RQ_QURS.csv'))
  rqq_csv(rqq_resumir_predicoes(val,cenario,config),file.path(pasta,'Resumo_Predicoes_RQ_QURS.csv'))
  if(nrow(descartadas)) {
    tab<-as.data.frame(table(descartadas$tipo,descartadas$erro),stringsAsFactors=FALSE)
    tab<-tab[tab$Freq>0,,drop=FALSE];names(tab)<-c('tipo','erro','frequencia')
  } else tab<-data.frame(tipo=character(),erro=character(),frequencia=integer())
  rqq_csv(tab,file.path(pasta,'Resumo_Erros.csv'))
  invisible(controle)
}

# Run one configuration, recording accepted and rejected attempts.
rqq_executar_configuracao <- function(config,cenario,base,pasta,B=10000L,
  semente_base=20260817L,controle=list(),checkpoint_cada=100L,
  maximo_tentativas=1000000L,maximo_falhas_consecutivas=2000L,
  retomar=TRUE,assinatura_codigo='',max_novas_tentativas=Inf,
  gerador=rRQ_QURS,estimador=estimar_RQ_QURS) {
  rqq_assert(B>=1 && B==floor(B) && B<=maximo_tentativas,'B invalido.')
  rqq_assert(checkpoint_cada>=1 && maximo_tentativas<=1000000L &&
    maximo_falhas_consecutivas>=1 && max_novas_tentativas>=0,'Controles de execucao invalidos.')
  rqq_assert(nrow(config)==1 && config$n==nrow(base),'Configuracao incompativel com desenho.')
  X<-rqq_matrizes_cenario(cenario,base);nm<-rqq_nomes(X)
  verdade<-setNames(rqq_verdade(cenario),nm)
  predverd<-plogis(drop(X%*%verdade))
  rqq_parametros(predverd,config$tau)
  dir.create(pasta,recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(pasta,'.em_execucao')
  rqq_assert(dir.create(lock,showWarnings=FALSE),paste('Pasta bloqueada por outra execucao:',lock))
  writeLines(c(paste('pid',Sys.getpid()),paste('host',Sys.info()[['nodename']]),as.character(Sys.time())),
             file.path(lock,'responsavel.txt'))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  arquivo<-file.path(pasta,'Checkpoint_RQ_QURS.rds')
  assinatura<-list(versao='RQ_QURS_gosolnp_v1',config=config,cenario=cenario,
    X=X,B=as.integer(B),semente_base=semente_base,
    controle=modifyList(rqq_controle_gosolnp(),controle),maximo_tentativas=maximo_tentativas,
    R=as.character(getRversion()),Rsolnp=as.character(utils::packageVersion('Rsolnp')),
    codigo=assinatura_codigo)
  if(file.exists(arquivo)) {
    rqq_assert(retomar,'Ja existe checkpoint: use retomar=TRUE ou outra pasta.')
    estado<-rqq_ler_checkpoint(arquivo)
    rqq_assert(identical(estado$assinatura,assinatura),
      'Checkpoint incompativel: configuracao, desenho, B, sementes, versoes ou codigo mudaram. Use nova pasta.')
  } else estado<-list(assinatura=assinatura,validos=list(),tentativas=list(),falhas_consecutivas=0L)
  salvar<-function() {
    rqq_salvar_rds(estado,arquivo)
    rqq_exportar_configuracao(estado,pasta,config,cenario,X,verdade)
  }
  novas<-0L
  # Use separate seed ranges for sampling, optimization, and design generation.
  # A stride of one million preserves seeds when subsetting or resuming.
  while(length(estado$validos)<B && novas<max_novas_tentativas) {
    t<-length(estado$tentativas)+1L
    if(t>maximo_tentativas || estado$falhas_consecutivas>=maximo_falhas_consecutivas) {
      salvar();stop('Limite de tentativas/falhas consecutivas atingido. Resultado INCOMPLETO salvo.',call.=FALSE)
    }
    sg<-rqq_semente(semente_base+(config$id_configuracao-1)*1000000+t)
    so<-rqq_semente(1000000000+semente_base+(config$id_configuracao-1)*1000000+t)
    start<-proc.time()[['elapsed']];etapa<-'geracao';avisos<-character()
    fit<-tryCatch(withCallingHandlers({
      RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(sg)
      y<-gerador(X,cenario$beta,config$tau)
      rqq_assert(length(y)==config$n && all(is.finite(y)) && all(y>0 & y<1),
        'Gerador retornou amostra nao finita ou fora de (0,1).')
      etapa<-'estimacao'
      f<-estimador(y,X,config$tau,rseed=so,controle=controle)
      rqq_assert(is.list(f) && length(f$coef)==length(verdade) &&
        all(is.finite(c(f$coef,f$se,f$lower,f$upper))) &&
        all(f$se>0) && all(f$upper>f$lower),'Ajuste incompleto ou inferencia invalida.')
      # Check prediction variances and intervals before accepting a replication.
      ref<-rqq_pontos_referencia(cenario);pr<-rqq_predizer(f,ref$X)
      rqq_assert(all(is.finite(as.matrix(pr))) && all(pr$ls_logit>pr$li_logit),
        'Predicoes ou intervalos invalidos.')
      f
    },warning=function(w) {avisos<<-c(avisos,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e) e)
    tempo<-proc.time()[['elapsed']]-start;ok<-!inherits(fit,'error')
    tipo<-if(ok) 'OK' else if(inherits(fit,'erro_rqqurs')) fit$tipo else paste0('erro_',etapa)
    mensagem<-if(ok) '' else conditionMessage(fit)
    avisos<-unique(c(avisos,if(ok) fit$avisos else character()))
    row<-data.frame(tentativa=t,semente_amostra=sg,semente_otimizacao=so,valida=ok,
      tipo=tipo,erro=mensagem,avisos=paste(avisos,collapse=' | '),tempo=tempo,
      convergencia=if(ok) fit$convergencia else if(!is.null(fit$convergencia)) fit$convergencia else NA_integer_,
      gradiente_max=if(ok) fit$gradiente_max else if(!is.null(fit$gradiente)) fit$gradiente else NA_real_,
      rcond=if(ok) fit$rcond else if(!is.null(fit$rcond)) fit$rcond else NA_real_,stringsAsFactors=FALSE)
    estado$tentativas[[t]]<-row
    if(ok) {
      fit$tentativa<-t;fit$semente_amostra<-sg;fit$tempo<-tempo;fit$avisos<-avisos
      estado$validos[[length(estado$validos)+1L]]<-fit;estado$falhas_consecutivas<-0L
    } else estado$falhas_consecutivas<-estado$falhas_consecutivas+1L
    novas<-novas+1L
    if(t%%checkpoint_cada==0L) {
      salvar();message(sprintf('%s tau=%.1f n=%d: %d/%d accepted; %d attempts; %d discarded.',
        config$cenario,config$tau,config$n,length(estado$validos),B,t,t-length(estado$validos)))
    }
  }
  salvar();invisible(estado)
}
# Merge configuration-level summaries into study-wide CSV files.
rqq_compilar_resultados <- function(raiz) {
  arquivos<-list.files(file.path(raiz,'CENARIOS'),recursive=TRUE,full.names=TRUE)
  destinos<-c(Controle_Replicas.csv='Controle_Replicas_Geral.csv',
    Resumo_Monte_Carlo_RQ_QURS.csv='Resumo_Monte_Carlo_RQ_QURS_Geral.csv',
    Resumo_Predicoes_RQ_QURS.csv='Resumo_Predicoes_RQ_QURS_Geral.csv')
  for(nm in names(destinos)) {
    fs<-arquivos[basename(arquivos)==nm]
    dfs<-lapply(fs,function(f) tryCatch(utils::read.csv(f,stringsAsFactors=FALSE),error=function(e) NULL))
    dfs<-Filter(function(x) !is.null(x)&&nrow(x)>0,dfs)
    if(length(dfs)) rqq_csv(do.call(rbind,dfs),file.path(raiz,destinos[[nm]]))
  };invisible(raiz)
}




