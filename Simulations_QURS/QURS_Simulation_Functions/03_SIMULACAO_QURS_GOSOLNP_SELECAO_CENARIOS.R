# QURS Monte Carlo driver: retain B replications with complete inference.
# Run with Rscript or source this file in RStudio. Keep the helper file alongside it.
# See README.md for the full run and a one-configuration smoke test.
modo_execucao <- getOption("QURS.mode", "final") # Values: teste, piloto, intermediario, final.
n_replicacoes_personalizado <- NULL
semente_base <- 20261016L
# Optimizer settings are part of the checkpoint signature.
gosolnp_n_restarts <- 5L
gosolnp_n_sim <- 30L
gosolnp_margem <- 1e-8
gosolnp_controle <- list(trace=0,outer.iter=400L,inner.iter=800L,tol=1e-8,delta=1e-7)
nivel_confianca <- 0.95
usar_paralelismo <- getOption("QURS.parallel", TRUE)
maximo_nucleos <- 4L
retomar_execucao <- TRUE
tamanho_lote <- 50L
maximo_tentativas_por_configuracao <- 1000000L
# Use one scenario selector at a time. NULL for both selects all scenarios.
# From C2 onward: cenario_inicial <- 'C2'; cenarios_selecionados <- NULL.
# C1 and C3 only: cenario_inicial <- NULL; cenarios_selecionados <- c('C1','C3').
cenario_inicial <- NULL
cenarios_selecionados <- getOption("QURS.scenarios", NULL)
taus_selecionados <- getOption("QURS.tau", NULL)
ns_selecionados <- getOption("QURS.n", NULL)
versao_estudo <- 'gosolnp_v1_validas'

obter_pasta_QURS <- function() {
  arqs <- vapply(sys.frames(), function(x)
    if (is.null(x$ofile)) NA_character_ else as.character(x$ofile)[1], character(1))
  arqs <- arqs[!is.na(arqs)]
  if (length(arqs)) return(dirname(normalizePath(tail(arqs,1),winslash='/',mustWork=TRUE)))
  argumento <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if (length(argumento)) return(dirname(normalizePath(sub('^--file=','',argumento[1]),mustWork=TRUE)))
  normalizePath(getwd(),winslash='/',mustWork=TRUE)
}
pasta_script <- obter_pasta_QURS()
arquivo_funcoes <- file.path(pasta_script,'01_FUNCOES_QURS_GOSOLNP.R')
if (!file.exists(arquivo_funcoes)) stop('File not found: ',arquivo_funcoes)
source(arquivo_funcoes, local=FALSE)
if(!requireNamespace('Rsolnp',quietly=TRUE))
  stop("Instale install.packages('Rsolnp') antes de executar a simulacao.")
versao_Rsolnp <- as.character(utils::packageVersion('Rsolnp'))
for(x in list(gosolnp_n_restarts,gosolnp_n_sim))
  if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<1 || x!=floor(x))
    stop('gosolnp controls must be positive integers.')
if(gosolnp_n_sim<gosolnp_n_restarts) stop('n.sim deve ser >= n.restarts.')
if(length(gosolnp_margem)!=1L || !is.finite(gosolnp_margem) || gosolnp_margem<=0 || gosolnp_margem>=.1)
  stop('Invalid numerical boundary margin.')
replicacoes_por_modo <- c(teste=2L,piloto=20L,intermediario=1000L,final=5000L)
if (!modo_execucao %in% names(replicacoes_por_modo)) stop('Invalid mode. Use teste, piloto, intermediario, or final.')
n_replicacoes <- if (is.null(n_replicacoes_personalizado))
  unname(replicacoes_por_modo[modo_execucao]) else n_replicacoes_personalizado
for (valor in list(n_replicacoes,tamanho_lote,maximo_nucleos,maximo_tentativas_por_configuracao))
  if (!is.numeric(valor) || length(valor)!=1L || !is.finite(valor) || valor<1 || valor!=floor(valor))
    stop('Replications, batch size, core count, and attempt limit must be positive integers.')
.validar_prob_QURS(nivel_confianca,'nivel_confianca')
if (!is.numeric(semente_base) || length(semente_base)!=1L || !is.finite(semente_base) ||
    semente_base<0 || semente_base!=floor(semente_base) || semente_base>=.Machine$integer.max)
  stop('Invalid base seed.')
if (maximo_tentativas_por_configuracao<n_replicacoes) stop('The attempt limit is smaller than B.')

# Build the complete factorial design before applying any filters.
# Keep configuration IDs fixed so filtered runs use the same seeds.
mu_tau_valores <- c(0.1,0.3,0.5,0.7,0.9)
tau_valores <- c(0.10,0.25,0.50,0.75,0.90)
n_valores <- c(50L,100L,200L,400L)
grade_completa <- do.call(rbind,lapply(seq_along(mu_tau_valores),function(i) {
  g <- expand.grid(n=n_valores,tau=tau_valores,KEEP.OUT.ATTRS=FALSE)
  g$mu_tau <- mu_tau_valores[i]; g$cenario <- paste0('C',i); g
}))
row.names(grade_completa) <- NULL
grade_completa$id_configuracao <- seq_len(nrow(grade_completa))
normalizar_cenarios_QURS <- function(x,disponiveis) {
  if(is.numeric(x)) {
    if(any(!is.finite(x)) || any(x!=floor(x))) stop('Numeric scenario IDs must be integers.')
    x <- paste0('C',x)
  } else if(is.character(x)) x <- toupper(trimws(x)) else
    stop('Use scenario labels such as C3 or integer IDs such as 3.')
  if(!length(x) || anyNA(x) || any(!x %in% disponiveis))
    stop('Unknown scenario. Available scenarios: ',paste(disponiveis,collapse=', '))
  unique(x)
}
selecionar_grade_QURS <- function(completa,cenario_inicial=NULL,cenarios_selecionados=NULL) {
  disponiveis <- unique(completa$cenario)
  if(!is.null(cenario_inicial) && !is.null(cenarios_selecionados))
    stop('Use only one scenario selector: cenario_inicial or cenarios_selecionados.')
  escolhidos <- disponiveis
  if(!is.null(cenario_inicial)) {
    inicial <- normalizar_cenarios_QURS(cenario_inicial,disponiveis)
    if(length(inicial)!=1L) stop('cenario_inicial must identify exactly one scenario.')
    escolhidos <- disponiveis[seq.int(match(inicial,disponiveis),length(disponiveis))]
  }
  if(!is.null(cenarios_selecionados))
    escolhidos <- normalizar_cenarios_QURS(cenarios_selecionados,disponiveis)
  # Preserve row order and IDs when selecting a subset.
  completa[completa$cenario %in% escolhidos,,drop=FALSE]
}

grade_estudo <- selecionar_grade_QURS(grade_completa,cenario_inicial,cenarios_selecionados)
for (filtro in list(list('tau',taus_selecionados),
                    list('n',ns_selecionados))) {
  if (!is.null(filtro[[2]])) {
    if (any(!filtro[[2]] %in% grade_completa[[filtro[[1]]]])) stop('Unknown filter value for ',filtro[[1]])
    grade_estudo <- grade_estudo[grade_estudo[[filtro[[1]]]] %in% filtro[[2]],,drop=FALSE]
  }
}
if (!nrow(grade_estudo)) stop('The selected grid is empty.')
cat('Selected scenarios: ',paste(unique(grade_estudo$cenario),collapse=', '),'\n')
if (maximo_tentativas_por_configuracao*nrow(grade_completa)>=.Machine$integer.max-1)
  stop('The attempt limit exceeds the range of unique seeds.')
output_root <- getOption('QURS.output_root', file.path(pasta_script, 'RESULTADOS_MONTE_CARLO'))
pasta_resultados <- file.path(output_root,
  paste0(modo_execucao,'_',versao_estudo,'_R',formatC(n_replicacoes,width=5,flag='0')))
pasta_replicacoes <- file.path(pasta_resultados,'REPLICACOES')
dir.create(pasta_replicacoes,recursive=TRUE,showWarnings=FALSE)
tag_QURS <- function(x) gsub('\\.','',formatC(x,format='f',digits=1))
arquivo_configuracao <- function(cfg) file.path(pasta_replicacoes,
  paste0(cfg$cenario,'_tau_',tag_QURS(cfg$tau),'_n_',formatC(cfg$n,width=3,flag='0'),'.rds'))

salvar_rds_QURS <- function(objeto,caminho) {
  tmp <- paste0(caminho,'.tmp'); bak <- paste0(caminho,'.bak')
  saveRDS(objeto,tmp)
  # Use atomic replacement where supported.
  if (suppressWarnings(file.rename(tmp,caminho))) return(invisible(caminho))
  # Keep a backup when the platform cannot replace an existing file by renaming.
  if (file.exists(caminho) && !file.copy(caminho,bak,overwrite=TRUE))
    stop('Could not preserve backup: ',caminho)
  if (!file.copy(tmp,caminho,overwrite=TRUE)) stop('Could not write: ',caminho)
  unlink(tmp)
  invisible(caminho)
}
ler_checkpoint_QURS <- function(caminho) {
  if (file.exists(caminho)) {
    x <- tryCatch(readRDS(caminho),error=function(e) e)
    if (!inherits(x,'error')) return(x)
  }
  backup <- paste0(caminho,'.bak')
  if (file.exists(backup)) {
    message('Restoring the previous checkpoint: ',backup)
    return(readRDS(backup))
  }
  stop('Missing or unreadable checkpoint: ',caminho)
}
# Derive deterministic seeds from the full-grid ID and attempt number.
semente_tentativa_QURS <- function(id,tentativa) as.integer(1+
  (as.double(semente_base)+(as.double(tentativa)-1)*nrow(grade_completa)+id) %%
  (.Machine$integer.max-1))
linha_vazia_QURS <- function(cfg,tentativa) data.frame(
  id_configuracao=cfg$id_configuracao,cenario=cfg$cenario,tau=cfg$tau,n=cfg$n,
  mu_tau_verdadeiro=cfg$mu_tau,alpha_verdadeiro=alpha_QURS(cfg$mu_tau,cfg$tau),
  replica=NA_integer_,tentativa=tentativa,
  semente=semente_tentativa_QURS(cfg$id_configuracao,tentativa),
  semente_otimizacao=.semente_otimizacao_QURS(semente_tentativa_QURS(cfg$id_configuracao,tentativa)),
  metodo_otimizacao='Rsolnp::gosolnp',versao_Rsolnp=versao_Rsolnp,
  n_restarts=gosolnp_n_restarts,n_sim=gosolnp_n_sim,
  mu_tau_estimado=NA_real_,alpha_estimado=NA_real_,se_mu_tau=NA_real_,se_alpha=NA_real_,
  logLik_verdadeiro=NA_real_,logLik_estimado=NA_real_,informacao_log_alpha=NA_real_,
  escore_log_alpha=NA_real_,
  ic_wald_mu_inferior=NA_real_,ic_wald_mu_superior=NA_real_,
  ic_logit_mu_inferior=NA_real_,ic_logit_mu_superior=NA_real_,
  cobertura_wald_mu_tau=NA,cobertura_logit_mu_tau=NA,
  convergencia=NA_integer_,hessiana_positiva=FALSE,inferencia_regular=FALSE,
  estimativa_valida=FALSE,inferencia_wald_valida=FALSE,inferencia_logit_valida=FALSE,
  aceita=FALSE,status='PENDENTE',etapa='',motivo_descarte='',mensagem='',avisos='',
  tempo_segundos=NA_real_,stringsAsFactors=FALSE)

# Generate, fit, and validate one attempt. Record warnings and failures so
# the driver can replace rejected attempts and report the discard rate.
executar_tentativa_QURS <- function(tentativa,cfg,nivel_confianca) {
  linha <- linha_vazia_QURS(cfg,tentativa)
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  set.seed(linha$semente)
  inicio <- proc.time()['elapsed']; avisos <- character(); etapa <- 'GERACAO'
  resultado <- tryCatch(withCallingHandlers({
    y <- rQURS(cfg$n,cfg$mu_tau,cfg$tau)
    if (length(y)!=cfg$n || any(!is.finite(y)) || any(y<=0 | y>=1))
      stop('Sample contains nonfinite values or values outside (0, 1).')
    etapa <- 'VEROSSIMILHANCA_VERDADEIRA'
    linha$logLik_verdadeiro <- logLik_QURS(y,cfg$mu_tau,cfg$tau)
    if (!is.finite(linha$logLik_verdadeiro)) stop('Nonfinite log-likelihood at the true parameter.')
    etapa <- 'ESTIMACAO'
    fit <- estimar_QURS(y,cfg$tau,controle=gosolnp_controle,
      rseed=linha$semente_otimizacao,n.restarts=gosolnp_n_restarts,
      n.sim=gosolnp_n_sim,margem_numerica=gosolnp_margem)
    linha$mu_tau_estimado <- unname(fit$estimativas['mu_tau'])
    linha$alpha_estimado <- fit$alpha_estimado
    linha$se_mu_tau <- unname(fit$erro_padrao['mu_tau']); linha$se_alpha <- fit$se_alpha
    linha$convergencia <- fit$convergencia; linha$hessiana_positiva <- fit$hessiana_positiva
    linha$inferencia_regular <- fit$inferencia_regular
    linha$logLik_estimado <- fit$logLik
    linha$informacao_log_alpha <- fit$informacao_log_alpha
    linha$escore_log_alpha <- fit$escore_log_alpha
    if (!isTRUE(fit$convergencia==0L)) stop('FALHA_CONVERGENCIA: ',fit$mensagem)
    mu <- linha$mu_tau_estimado
    if (!is.finite(mu) || mu<=0 || mu>=1 || !is.finite(fit$alpha_estimado) || fit$alpha_estimado<=0)
      stop('ESTIMATIVA_INVALIDA_OU_FRONTEIRA_NUMERICA')
    linha$estimativa_valida <- TRUE
    if (!is.finite(fit$logLik) || fit$logLik<linha$logLik_verdadeiro-1e-6)
      stop('MAXIMIZACAO_INVALIDA: fitted log-likelihood is below the value at the true parameter.')
    etapa <- 'HESSIANA'
    if (!isTRUE(fit$hessiana_positiva)) stop('HESSIANA_INVALIDA')
    if (!isTRUE(fit$inferencia_regular) || !is.finite(linha$se_mu_tau) || linha$se_mu_tau<=0 ||
        !is.finite(linha$se_alpha) || linha$se_alpha<=0) stop('ERRO_PADRAO_INVALIDO')
    etapa <- 'INTERVALOS'
    w <- confint_QURS(fit,nivel_confianca,'wald')
    l <- confint_QURS(fit,nivel_confianca,'logit')
    valido <- function(x) is.matrix(x) && identical(dim(x),c(1L,2L)) &&
      all(is.finite(x)) && x[1,1]<x[1,2]
    linha$inferencia_wald_valida <- valido(w)
    linha$inferencia_logit_valida <- valido(l)
    if (!valido(w) || !valido(l)) stop('INTERVALO_INVALIDO')
    linha$ic_wald_mu_inferior <- w[1,1]; linha$ic_wald_mu_superior <- w[1,2]
    linha$ic_logit_mu_inferior <- l[1,1]; linha$ic_logit_mu_superior <- l[1,2]
    linha$cobertura_wald_mu_tau <- w[1,1]<=cfg$mu_tau && cfg$mu_tau<=w[1,2]
    linha$cobertura_logit_mu_tau <- l[1,1]<=cfg$mu_tau && cfg$mu_tau<=l[1,2]
    linha$aceita <- TRUE; linha$status <- 'OK'; linha
  },warning=function(w) {
    avisos <<- c(avisos,paste0(etapa,': ',conditionMessage(w)))
    invokeRestart('muffleWarning')
  }),error=function(e) {
    linha$status <- paste0('ERRO_',etapa)
    linha$motivo_descarte <- linha$status
    linha$mensagem <- conditionMessage(e)
    linha$aceita <- FALSE
    linha
  })
  if (length(avisos)) {
    resultado$aceita <- FALSE
    resultado$motivo_descarte <- paste(c(resultado$motivo_descarte,'AVISO_R')[
      nzchar(c(resultado$motivo_descarte,'AVISO_R'))],collapse=' | ')
    if (resultado$status=='OK') resultado$status <- 'AVISO_R'
  }
  resultado$etapa <- etapa; resultado$avisos <- paste(unique(avisos),collapse=' | ')
  resultado$tempo_segundos <- unname(proc.time()['elapsed']-inicio)
  resultado
}

resumo_QURS <- function(dados,descartadas,cfg,tentativas) {
  v <- dados$mu_tau_estimado; se <- dados$se_mu_tau; m <- mean(v); dp <- sd(v)
  data.frame(id_configuracao=cfg$id_configuracao,cenario=cfg$cenario,tau=cfg$tau,n=cfg$n,
    mu_tau_verdadeiro=cfg$mu_tau,nivel_confianca=nivel_confianca,replicas_solicitadas=n_replicacoes,replicas_salvas=nrow(dados),
    estimativas_validas=nrow(dados),intervalos_wald_validos=sum(dados$inferencia_wald_valida),
    intervalos_logit_validos=sum(dados$inferencia_logit_valida),
    tentativas_totais=tentativas,replicas_descartadas=NROW(descartadas),
    taxa_descarte=NROW(descartadas)/tentativas,
    mu_media=m,mu_vies=m-cfg$mu_tau,mu_vies_relativo_percentual=100*(m-cfg$mu_tau)/cfg$mu_tau,
    mu_desvio_padrao=dp,mu_erro_padrao_medio=mean(se),
    mu_eqm=mean((v-cfg$mu_tau)^2),mu_raiz_eqm=sqrt(mean((v-cfg$mu_tau)^2)),
    mu_assimetria=if(length(v)>=3 && is.finite(dp) && dp>0) mean(((v-m)/dp)^3) else NA_real_,
    mu_curtose=if(length(v)>=4 && is.finite(dp) && dp>0) mean(((v-m)/dp)^4) else NA_real_,
    mu_cobertura_wald=mean(dados$cobertura_wald_mu_tau),
    mu_nao_cobertura_wald=1-mean(dados$cobertura_wald_mu_tau),
    mu_cobertura_logit=mean(dados$cobertura_logit_mu_tau),
    mu_nao_cobertura_logit=1-mean(dados$cobertura_logit_mu_tau),
    mu_amplitude_media_wald=mean(dados$ic_wald_mu_superior-dados$ic_wald_mu_inferior),
    mu_amplitude_media_logit=mean(dados$ic_logit_mu_superior-dados$ic_logit_mu_inferior),
    tempo_medio_por_replica=mean(dados$tempo_segundos),
    tempo_total_das_replicas=sum(dados$tempo_segundos),
    tempo_total_tentativas=sum(dados$tempo_segundos)+sum(descartadas$tempo_segundos))
}

assinatura_configuracao_QURS <- function(cfg) list(versao=versao_estudo,cfg=cfg,
 B=n_replicacoes,semente=semente_base,grade_sementes=grade_completa,nivel=nivel_confianca,
 funcoes_md5=unname(tools::md5sum(arquivo_funcoes)),criterio='completa_sem_avisos_v1',
 otimizacao=list(metodo='Rsolnp::gosolnp',versao_Rsolnp=versao_Rsolnp,
  n_restarts=gosolnp_n_restarts,n_sim=gosolnp_n_sim,margem=gosolnp_margem,
  controle=gosolnp_controle))
controle_estado_QURS <- function(cfg,estado) data.frame(cfg,
 replicas_solicitadas=n_replicacoes,replicas_validas=NROW(estado$validas),
 replicas_descartadas=NROW(estado$descartadas),tentativas=estado$tentativas,
 taxa_descarte=if(estado$tentativas) NROW(estado$descartadas)/estado$tentativas else 0,
 completo=NROW(estado$validas)==n_replicacoes)
carregar_acumulados_QURS <- function() {
  vazio <- list(resumo=NULL,controle=NULL)
  f <- file.path(pasta_resultados,'Resumo_Monte_Carlo_QURS.csv')
  if(!isTRUE(retomar_execucao) || !file.exists(f)) return(vazio)
  anterior <- read.csv(f,stringsAsFactors=FALSE)
  manter <- logical(nrow(anterior)); controles <- list()
  for(k in seq_len(nrow(anterior))) {
    j <- match(anterior$id_configuracao[k],grade_completa$id_configuracao)
    if(is.na(j)) next
    cfg <- grade_completa[j,,drop=FALSE]
    ck <- file.path(pasta_resultados,'CENARIOS',cfg$cenario,
      paste0('tau_',tag_QURS(cfg$tau),'_n_',cfg$n),'Checkpoint.rds')
    if(!file.exists(ck) && !file.exists(paste0(ck,'.bak'))) next
    e <- ler_checkpoint_QURS(ck)
    if(!identical(e$assinatura,assinatura_configuracao_QURS(cfg))) next
    manter[k] <- NROW(e$validas)==n_replicacoes && all(e$validas$aceita)
    if(manter[k]) controles[[length(controles)+1L]] <- controle_estado_QURS(cfg,e)
  }
  list(resumo=anterior[manter,,drop=FALSE],controle=if(length(controles)) do.call(rbind,controles) else NULL)
}
unir_resultados_QURS <- function(anterior,novo) {
  if(!is.null(anterior) && nrow(anterior))
    novo <- rbind(anterior[!anterior$id_configuracao %in% novo$id_configuracao,,drop=FALSE],novo)
  novo <- novo[order(novo$id_configuracao),,drop=FALSE];row.names(novo) <- NULL
  novo
}

# Process selected configurations in batches and checkpoint each batch.
executar_estudo_QURS <- function() {
  nucleos <- 1L
  if (isTRUE(usar_paralelismo)) {
    nc <- parallel::detectCores(); if (!is.finite(nc)) nc <- 1L
    nucleos <- max(1L,min(maximo_nucleos,nc-1L))
  }
  cluster <- NULL
  if (nucleos>1L) {
    cluster <- parallel::makeCluster(nucleos)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
    parallel::clusterExport(cluster,c('arquivo_funcoes','semente_base','grade_completa','versao_Rsolnp',
      'gosolnp_n_restarts','gosolnp_n_sim','gosolnp_margem','gosolnp_controle'),
                            envir=environment())
    invisible(parallel::clusterEvalQ(cluster,{
      source(arquivo_funcoes,local=FALSE)
      if(!requireNamespace('Rsolnp',quietly=TRUE)) stop('Rsolnp is unavailable on the parallel worker.')
      NULL
    }))
    parallel::clusterExport(cluster,c('semente_tentativa_QURS','linha_vazia_QURS',
      'executar_tentativa_QURS'),envir=environment())
  }
  inicio <- proc.time()['elapsed']; resumos <- list(); controles <- list()
  acumulados <- carregar_acumulados_QURS()
  ids_planejados <- unique(c(acumulados$resumo$id_configuracao,grade_estudo$id_configuracao))
  planejada <- grade_completa[grade_completa$id_configuracao %in% ids_planejados,,drop=FALSE]
  write.csv(planejada,file.path(pasta_resultados,'Grade_Estudo_QURS.csv'),row.names=FALSE)
  write.csv(grade_estudo,file.path(pasta_resultados,'Grade_Selecao_Atual_QURS.csv'),row.names=FALSE)
  write.csv(grade_completa,file.path(pasta_resultados,'Grade_Completa_QURS.csv'),row.names=FALSE)
  capture.output(list(versao=versao_estudo,B=n_replicacoes,semente=semente_base,
    nivel=nivel_confianca,nucleos=nucleos,criterio='Complete estimation and inference, with no warnings',
    cenario_inicial=cenario_inicial,cenarios_selecionados=cenarios_selecionados,
    taus_selecionados=taus_selecionados,ns_selecionados=ns_selecionados,
    otimizacao=list(metodo='Rsolnp::gosolnp',versao=versao_Rsolnp,
      n_restarts=gosolnp_n_restarts,n_sim=gosolnp_n_sim,margem=gosolnp_margem,
      controle=gosolnp_controle),
    sessionInfo=sessionInfo()),file=file.path(pasta_resultados,'Configuracao_Execucao.txt'))
  for (i in seq_len(nrow(grade_estudo))) {
    cfg <- grade_estudo[i,,drop=FALSE]
    pasta <- file.path(pasta_resultados,'CENARIOS',cfg$cenario,
      paste0('tau_',tag_QURS(cfg$tau),'_n_',cfg$n))
    dir.create(pasta,recursive=TRUE,showWarnings=FALSE)
    checkpoint <- file.path(pasta,'Checkpoint.rds')
    assinatura <- assinatura_configuracao_QURS(cfg)
    estado <- list(assinatura=assinatura,validas=NULL,descartadas=NULL,tentativas=0L)
    if (isTRUE(retomar_execucao) && (file.exists(checkpoint)||file.exists(paste0(checkpoint,'.bak')))) {
      estado <- ler_checkpoint_QURS(checkpoint)
      if (!identical(estado$assinatura,assinatura)) stop('Incompatible checkpoint: ',checkpoint)
    }
    controle <- NULL
    salvar <- function() {
      salvar_rds_QURS(estado,checkpoint) # Use the checkpoint as the source for resuming the run.
      vazia <- linha_vazia_QURS(cfg,1L)[0,,drop=FALSE]
      v <- if(is.null(estado$validas)) vazia else estado$validas
      d <- if(is.null(estado$descartadas)) vazia else estado$descartadas
      for (par in list(list(v,'Replicas_Validas'),list(d,'Replicas_Descartadas'))) {
        salvar_rds_QURS(par[[1]],file.path(pasta,paste0(par[[2]],'.rds')))
        write.csv(par[[1]],file.path(pasta,paste0(par[[2]],'.csv')),row.names=FALSE,fileEncoding='UTF-8')
      }
      salvar_rds_QURS(v,arquivo_configuracao(cfg))
      tab <- data.frame(status=character(),motivo_descarte=character(),quantidade=integer())
      if (nrow(d)) tab <- aggregate(list(quantidade=rep(1L,nrow(d))),
        d[c('status','motivo_descarte')],sum)
      write.csv(tab,file.path(pasta,'Contagem_Motivos_Descarte.csv'),row.names=FALSE)
      controle <<- controle_estado_QURS(cfg,estado)
      write.csv(controle,file.path(pasta,'Controle_Replicas.csv'),row.names=FALSE)
    }
    cat(sprintf('[%d/%d] %s: mu=%s tau=%s n=%s\n',i,nrow(grade_estudo),
                cfg$cenario,cfg$mu_tau,cfg$tau,cfg$n))
    repeat {
      nv <- NROW(estado$validas)
      if (nv==n_replicacoes) break
      if (estado$tentativas>=maximo_tentativas_por_configuracao) {
        salvar()
        stop('Attempt limit reached. Study incomplete: ',pasta,
             '. Review discarded attempts before increasing the limit and resuming.')
      }
      quantidade <- min(tamanho_lote,n_replicacoes-nv,
                        maximo_tentativas_por_configuracao-estado$tentativas)
      ids <- estado$tentativas+seq_len(quantidade)
      linhas <- if(nucleos>1L) parallel::parLapplyLB(cluster,ids,executar_tentativa_QURS,
        cfg=cfg,nivel_confianca=nivel_confianca) else
        lapply(ids,executar_tentativa_QURS,cfg=cfg,nivel_confianca=nivel_confianca)
      linhas <- do.call(rbind,linhas)
      boas <- linhas[linhas$aceita,,drop=FALSE]; ruins <- linhas[!linhas$aceita,,drop=FALSE]
      if (nrow(boas)) {
        boas$replica <- nv+seq_len(nrow(boas)); estado$validas <- rbind(estado$validas,boas)
      }
      if (nrow(ruins)) estado$descartadas <- rbind(estado$descartadas,ruins)
      estado$tentativas <- max(ids)
      salvar()
      cat(sprintf('  accepted=%d/%d, discarded=%d, attempts=%d\n',
        NROW(estado$validas),n_replicacoes,NROW(estado$descartadas),estado$tentativas))
    }
    salvar()
    stopifnot(NROW(estado$validas)==n_replicacoes,all(estado$validas$aceita),
      !anyDuplicated(estado$validas$tentativa),
      estado$tentativas==NROW(estado$validas)+NROW(estado$descartadas))
    resumo <- resumo_QURS(estado$validas,estado$descartadas,cfg,estado$tentativas)
    write.csv(resumo,file.path(pasta,'Resumo_Monte_Carlo_QURS.csv'),row.names=FALSE)
    resumos[[i]] <- resumo; controles[[i]] <- controle
    resumo_selecao <- do.call(rbind,resumos); controle_selecao <- do.call(rbind,controles)
    write.csv(resumo_selecao,file.path(pasta_resultados,'Resumo_Selecao_Atual_QURS.csv'),row.names=FALSE)
    write.csv(controle_selecao,file.path(pasta_resultados,'Controle_Selecao_Atual_QURS.csv'),row.names=FALSE)
    write.csv(unir_resultados_QURS(acumulados$resumo,resumo_selecao),file.path(pasta_resultados,'Resumo_Monte_Carlo_QURS.csv'),row.names=FALSE)
    write.csv(unir_resultados_QURS(acumulados$controle,controle_selecao),file.path(pasta_resultados,'Controle_Global_QURS.csv'),row.names=FALSE)
  }
  resumo <- do.call(rbind,resumos)
  diagnostico <- data.frame(configuracoes=nrow(resumo),replicas_por_configuracao=n_replicacoes,
    replicas_validas=sum(resumo$replicas_salvas),tentativas_totais=sum(resumo$tentativas_totais),
    descartadas=sum(resumo$replicas_descartadas),nucleos=nucleos,
    tempo_total_observado_segundos=unname(proc.time()['elapsed']-inicio),
    tempo_acumulado_tentativas_segundos=sum(resumo$tempo_total_tentativas),
    estimativa_final_serial_horas=sum(resumo$tempo_total_tentativas)/sum(resumo$replicas_salvas)*
      nrow(grade_estudo)*5000/3600)
  diagnostico$estimativa_final_paralela_ideal_horas <- diagnostico$estimativa_final_serial_horas/nucleos
  write.csv(diagnostico,file.path(pasta_resultados,'Diagnostico_Tempo_QURS.csv'),row.names=FALSE)
  list(resumo=resumo,diagnostico_tempo=diagnostico,pasta_resultados=pasta_resultados)
}

# To load definitions without running the study: options(QURS.executar=FALSE).
if (isTRUE(getOption('QURS.executar',TRUE))) {
  resultado_estudo <- executar_estudo_QURS()
  cat('\nStudy complete: B accepted replications in each selected configuration.\n')
  cat('Results: ',normalizePath(pasta_resultados,winslash='/'), '\n')
}


