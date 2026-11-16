# Source in RStudio or run with Rscript. Requires Rsolnp.
# Run script 04 for numerical checks before the full study.

# ----------------------- EDITABLE SETTINGS ------------------------------
modo <- getOption('QRQURS.mode', 'final')                        # Modes: teste (20), piloto (200), final (5000).
B_por_modo <- c(teste=20L,piloto=200L,final=5000L)
taus <- c(.1,.25,.5,.75,.9)
tamanhos_amostrais <- c(50L,100L,200L,400L)
cenario_inicial <- NULL                 # Example: 'R2' selects R2 and later scenarios.
cenarios_selecionados <- getOption('QRQURS.scenarios', NULL)           # Example: c('R1','R3').
incluir_cenario_adicional <- TRUE       # R4 uses correlated covariates.
semente_base <- 20260817L
semente_desenho <- 220260817L
numero_nucleos <- getOption('QRQURS.workers', 1L)                    # 1 = serial; >1 = PSOCK workers across configurations.
checkpoint_cada <- 100L
maximo_tentativas_por_configuracao <- 1000000L
maximo_falhas_consecutivas <- 2000L
retomar <- TRUE
# Use a different suffix when changing the study design, settings, or software.
sufixo_experimento <- getOption('QRQURS.suffix', '')                # Example: '_study2'.
controle_otimizacao <- list(n.restarts=5L,n.sim=200L,sd_inicial=1.5,
  limite_beta=12,expansoes=2L,tol_gradiente=1e-4,tol_rcond=1e-10,
  margem_preditor=1e-10,control=list(trace=0,outer.iter=400L,inner.iter=800L,tol=1e-8,delta=1e-7))
# -----------------------------------------------------------------------

.rqq_pasta <- local({
  args<-commandArgs(trailingOnly=FALSE);ff<-grep('^--file=',args,value=TRUE)
  of<-Filter(function(x) !is.null(x),lapply(sys.frames(),function(x) x$ofile))
  f<-if(length(of)) tail(of,1)[[1]] else if(length(ff)) sub('^--file=','',ff[1]) else NULL
  if(is.null(f)) {
    if(file.exists('01_FUNCOES_RQ_QURS_GOSOLNP.R')) getwd() else
      stop('Run the entire file with source() or Rscript.')
  } else dirname(normalizePath(f,mustWork=TRUE))
})
source(file.path(.rqq_pasta,'01_FUNCOES_RQ_QURS_GOSOLNP.R'),encoding='UTF-8')
source(file.path(.rqq_pasta,'02_MOTOR_MONTE_CARLO_RQ_QURS.R'),encoding='UTF-8')
if(!requireNamespace('Rsolnp',quietly=TRUE)) stop("Install Rsolnp with install.packages('Rsolnp').")
rqq_assert(modo%in%names(B_por_modo),'Invalid mode.');B<-unname(B_por_modo[[modo]])
rqq_assert(numero_nucleos>=1 && numero_nucleos==floor(numero_nucleos),'Invalid worker count.')
rqq_assert(!grepl('[/\\\\]',sufixo_experimento),'Sufixo nao pode conter separadores de caminho.')

# R1 includes an intercept; R2 and R3 omit it.
# See rqq_cenarios() in script 02 for the coefficients and design.
cenarios <- rqq_cenarios(incluir_cenario_adicional)
# Example: cenarios$R1$beta <- c(.5,.3,.2)

grade_completa<-rqq_grade(cenarios,taus,tamanhos_amostrais)
grade_estudo<-rqq_selecionar_grade(grade_completa,cenario_inicial,cenarios_selecionados)
rqq_assert(semente_base+nrow(grade_completa)*1000000+1000000000<.Machine$integer.max,
  'Grade excedeu a faixa reservada de sementes; reveja o desenho do estudo.')
desenhos<-rqq_desenhos_fixos(tamanhos_amostrais,semente_desenho)
raiz_resultados<-file.path(.rqq_pasta,'RESULTADOS_MONTE_CARLO',
  paste0(modo,'_rq_qurs_gosolnp_v1_validas_R',B,sufixo_experimento))
dir.create(raiz_resultados,recursive=TRUE,showWarnings=FALSE)
catalogo<-rqq_linhas(lapply(cenarios,function(c) data.frame(cenario=c$cenario,tipo=c$tipo,
  intercepto=c$intercepto,beta=paste(c$beta,collapse=';'),rho_latente=c$rho,
  descricao=c$descricao,stringsAsFactors=FALSE)))
hash_codigo<-unname(tools::md5sum(file.path(.rqq_pasta,c('01_FUNCOES_RQ_QURS_GOSOLNP.R',
                                            '02_MOTOR_MONTE_CARLO_RQ_QURS.R'))))
meta<-list(versao='RQ_QURS_gosolnp_v1',cenarios=cenarios,grade_completa=grade_completa,desenhos=desenhos,
  semente_desenho=semente_desenho,semente_base=semente_base,controle=controle_otimizacao,
  B=B,modo=modo,hash_codigo=hash_codigo,R=as.character(getRversion()),
  Rsolnp=as.character(utils::packageVersion('Rsolnp')))
arquivo_meta<-file.path(raiz_resultados,'Metadados_Estudo.rds')
if(file.exists(arquivo_meta)) rqq_assert(identical(readRDS(arquivo_meta),meta),
  'Pasta contem outro estudo. Mude sufixo_experimento ou restaure os controles originais.') else saveRDS(meta,arquivo_meta)
rqq_csv(grade_completa,file.path(raiz_resultados,'Grade_Completa.csv'))
rqq_csv(grade_estudo,file.path(raiz_resultados,'Grade_Selecionada_Ultima_Execucao.csv'))
rqq_csv(catalogo,file.path(raiz_resultados,'Cenarios.csv'))
capture.output(sessionInfo(),file=file.path(raiz_resultados,'SessionInfo.txt'))
for(c in cenarios) for(n in tamanhos_amostrais) {
  pastaX<-file.path(raiz_resultados,'DESENHOS',c$cenario)
  dir.create(pastaX,recursive=TRUE,showWarnings=FALSE)
  X<-rqq_matrizes_cenario(c,desenhos[[as.character(n)]])
  rqq_csv(as.data.frame(X),file.path(pastaX,paste0('X_n_',n,'.csv')))
}
rqq_tarefa <- function(i) {
  cfg<-grade_estudo[i,,drop=FALSE];c<-cenarios[[cfg$cenario]]
  pasta<-file.path(raiz_resultados,'CENARIOS',cfg$cenario,
    paste0('tau_',rqq_tag_tau(cfg$tau),'_n_',cfg$n))
  tryCatch({
    estado<-rqq_executar_configuracao(cfg,c,desenhos[[as.character(cfg$n)]],pasta,B,
      semente_base,controle_otimizacao,checkpoint_cada,maximo_tentativas_por_configuracao,
      maximo_falhas_consecutivas,retomar,hash_codigo)
    data.frame(cfg,completo=length(estado$validos)==B,erro='',stringsAsFactors=FALSE)
  },error=function(e) data.frame(cfg,completo=FALSE,erro=conditionMessage(e),stringsAsFactors=FALSE))
}
message(sprintf('QR-QURS gosolnp: %d configurations, %d valid replications per configuration.',nrow(grade_estudo),B))
if(numero_nucleos==1L) {
  resultados_execucao<-lapply(seq_len(nrow(grade_estudo)),rqq_tarefa)
} else {
  cl<-parallel::makePSOCKcluster(min(numero_nucleos,nrow(grade_estudo)),outfile='')
  resultados_execucao<-tryCatch({
    parallel::clusterExport(cl,'.rqq_pasta',envir=environment())
    parallel::clusterEvalQ(cl,{
      source(file.path(.rqq_pasta,'01_FUNCOES_RQ_QURS_GOSOLNP.R'),encoding='UTF-8')
      source(file.path(.rqq_pasta,'02_MOTOR_MONTE_CARLO_RQ_QURS.R'),encoding='UTF-8')
    })
    parallel::clusterExport(cl,c('grade_estudo','cenarios','desenhos','raiz_resultados','B',
      'semente_base','controle_otimizacao','checkpoint_cada','maximo_tentativas_por_configuracao',
      'maximo_falhas_consecutivas','retomar','hash_codigo','rqq_tarefa'),envir=environment())
    parallel::parLapply(cl,seq_len(nrow(grade_estudo)),rqq_tarefa)
  },finally=parallel::stopCluster(cl))
}
status<-rqq_linhas(resultados_execucao)
rqq_csv(status,file.path(raiz_resultados,'Status_Ultima_Execucao.csv'))
rqq_compilar_resultados(raiz_resultados)
if(any(!status$completo)) stop('Some configurations are incomplete. Check Status_Ultima_Execucao.csv before resuming.')
message('All selected configurations reached the requested replication count. Results: ',raiz_resultados)




