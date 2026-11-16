# Generate CSV and LaTeX tables, boxplots, and sampling-density plots.
# Read completed simulation outputs without refitting any models.
# An optional command-line argument selects a completed results directory.
args <- commandArgs(trailingOnly = TRUE)
pasta_resultados_manual <- if (length(args)) args[1L] else NULL
# Center each configuration on its empirical mean, as in the article.
# Other options: 'parametro_verdadeiro' or 'erro_padrao'.
centragem_padronizada <- 'media_empirica'
limites_z <- c(-4,4)
mostrar_legenda_no_pdf <- FALSE # Full captions are also written to a text file.

arqs <- vapply(sys.frames(),function(x) if(is.null(x$ofile)) NA_character_ else x$ofile[1],character(1))
arqs <- arqs[!is.na(arqs)]
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
pasta_script <- if (length(arqs)) {
  dirname(normalizePath(tail(arqs, 1)))
} else if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else getwd()
pasta_resultados <- if(is.null(pasta_resultados_manual)) file.path(pasta_script,
  'RESULTADOS_MONTE_CARLO','final_gosolnp_v1_validas_R05000') else pasta_resultados_manual
if (!centragem_padronizada %in% c('parametro_verdadeiro','media_empirica','erro_padrao'))
  stop('Unknown standardization method.')
arquivo_resumo <- file.path(pasta_resultados,'Resumo_Monte_Carlo_QURS.csv')
if (!file.exists(arquivo_resumo)) stop('Summary file not found: ',arquivo_resumo)
resumo <- read.csv(arquivo_resumo,stringsAsFactors=FALSE)
grade <- read.csv(file.path(pasta_resultados,'Grade_Estudo_QURS.csv'),stringsAsFactors=FALSE)
if (nrow(resumo)!=nrow(grade) || !setequal(resumo$id_configuracao,grade$id_configuracao))
  stop('Complete the selected simulation before generating plots.')
if (any(resumo$replicas_salvas!=resumo$replicas_solicitadas)) stop('Inconsistent replication counts.')
saida <- file.path(pasta_resultados,'TABELAS_GRAFICOS')
for (p in c('TABELAS','BOXPLOTS','DENSIDADES_ORIGINAIS','DENSIDADES_PADRONIZADAS','METRICAS'))
  dir.create(file.path(saida,p),recursive=TRUE,showWarnings=FALSE)
fmt <- function(x) gsub('\\.','',formatC(x,format='f',digits=1))
cores <- c('#24649B','#BD4F28','#32805B','#825BA0')
comb <- unique(resumo[c('cenario','tau')]); diagnosticos <- list(); legendas <- character()
num <- function(x) ifelse(is.finite(x),formatC(x,format='f',digits=5),'--')
# LaTeX fragments require no extra packages. Use \input{} in the main document.
for(i in seq_len(nrow(comb))) {
  sub <- resumo[resumo$cenario==comb$cenario[i] & resumo$tau==comb$tau[i],]
  sub <- sub[order(sub$n),]; tag <- paste0(comb$cenario[i],'_tau_',fmt(comb$tau[i]))
  write.csv(sub,file.path(saida,'TABELAS',paste0('Resumo_',tag,'.csv')),row.names=FALSE)
  tex <- c('% Results are conditional on replication acceptance; see the discard rate.',
    '\\begin{tabular}{rrrrrrrr}', '\\hline',
    '$n$ & Mean & Bias & SD & Mean SE & MSE & CP Wald & CP logit \\\\', '\\hline')
  for(j in seq_len(nrow(sub))) tex <- c(tex,paste0(paste(c(sub$n[j],
    num(unlist(sub[j,c('mu_media','mu_vies','mu_desvio_padrao','mu_erro_padrao_medio',
                       'mu_eqm','mu_cobertura_wald','mu_cobertura_logit')]))),collapse=' & '),' \\\\'))
  writeLines(c(tex,'\\hline','\\end{tabular}'),file.path(saida,'TABELAS',paste0('Tabela_',tag,'.tex')))
  amostras <- lapply(seq_len(nrow(sub)),function(j) {
    f <- file.path(pasta_resultados,'REPLICACOES',paste0(sub$cenario[j],
      '_tau_',fmt(sub$tau[j]),'_n_',formatC(sub$n[j],width=3,flag='0'),'.rds'))
    d <- readRDS(f)
    if (nrow(d)!=sub$replicas_solicitadas[j] || !all(d$aceita) ||
        any(!is.finite(d$mu_tau_estimado)) || anyDuplicated(d$tentativa))
      stop('Invalid replications in ',f)
    d
  })
  valores <- lapply(amostras,function(d) d$mu_tau_estimado)
  names(valores) <- paste0('n = ',sub$n)
  pdf(file.path(saida,'BOXPLOTS',paste0('Boxplot_mu_',tag,'.pdf')),width=6.4,height=4.7)
  par(mar=c(4.1,4.3,.7,.7),las=1)
  boxplot(valores,col=adjustcolor(cores,alpha.f=.55),xlab='',ylab=expression(hat(mu)[tau]),
          outline=TRUE,cex=.5)
  abline(h=sub$mu_tau_verdadeiro[1],col='#333333',lty=2,lwd=1.4)
  dev.off()
  padronizadas <- lapply(seq_along(amostras),function(j) {
    v <- valores[[j]]
    if (length(v)<2 || !is.finite(sd(v)) || sd(v)<=0)
      stop('Density estimation requires at least two varying estimates: ',tag)
    if (centragem_padronizada=='erro_padrao') return((v-sub$mu_tau_verdadeiro[j])/amostras[[j]]$se_mu_tau)
    centro <- if(centragem_padronizada=='media_empirica') mean(v) else sub$mu_tau_verdadeiro[j]
    (v-centro)/sd(v)
  })
  for(escala in c('original','padronizada')) {
    lista <- if(escala=='original') valores else padronizadas
    curvas <- lapply(lista,density,n=1024L)
    xlim <- if(escala=='original') pmax(0,pmin(1,range(unlist(lapply(curvas,function(d)d$x))))) else limites_z
    ylim <- range(0,unlist(lapply(curvas,function(d)d$y)),if(escala=='padronizada') dnorm(0) else 0)
    pasta <- if(escala=='original') 'DENSIDADES_ORIGINAIS' else 'DENSIDADES_PADRONIZADAS'
    pdf(file.path(saida,pasta,paste0('Densidade_mu_',tag,'.pdf')),width=6.4,height=4.7)
    par(mar=c(4.1,4.3,.7,.7),las=1)
    plot(NA,xlim=xlim,ylim=ylim,xlab=if(escala=='original') expression(hat(mu)[tau]) else 'z',
      ylab='Density')
    for(j in seq_along(curvas)) lines(curvas[[j]],col=cores[j],lwd=1.8,lty=j)
    if(escala=='original') abline(v=sub$mu_tau_verdadeiro[1],lty=2,col='#333333') else
      curve(dnorm(x),add=TRUE,lwd=1.6,lty=5,col='black',n=501)
    if(mostrar_legenda_no_pdf) legend('topright',legend=c(names(valores),
      if(escala=='padronizada') 'N(0,1)' else NULL),col=c(cores[seq_along(valores)],
      if(escala=='padronizada') 'black' else NULL),lty=c(seq_along(valores), if(escala=='padronizada') 5L else NULL),lwd=1.5,bty='n',cex=.8)
    dev.off()
  }
  for(j in seq_along(amostras)) {
    z <- padronizadas[[j]]
    diagnosticos[[length(diagnosticos)+1L]] <- data.frame(cenario=sub$cenario[j],tau=sub$tau[j],
      n=sub$n[j],centragem=centragem_padronizada,B=length(z),media_z=mean(z),dp_z=sd(z),
      assimetria_z=if(length(z)>=3) mean(((z-mean(z))/sd(z))^3) else NA_real_,
      curtose_z=if(length(z)>=4) mean(((z-mean(z))/sd(z))^4) else NA_real_,
      proporcao_fora_faixa=mean(z<limites_z[1]|z>limites_z[2]))
  }
  pdf(file.path(saida,'METRICAS',paste0('Metricas_',tag,'.pdf')),width=8,height=6)
  par(mfrow=c(2,2),mar=c(3.8,4.2,1.8,.8),las=1)
  plot(sub$n,sub$mu_vies,type='b',pch=16,col=cores[1],xlab='n',ylab='Bias');abline(h=0,lty=2)
  plot(sub$n,sub$mu_eqm,type='b',pch=16,col=cores[1],xlab='n',ylab='MSE')
  matplot(sub$n,sub[c('mu_cobertura_wald','mu_cobertura_logit')],type='b',pch=c(16,17),
    col=cores[1:2],lty=1:2,ylim=c(0,1),xlab='n',ylab='Coverage')
  abline(h=sub$nivel_confianca[1],lty=3);legend('bottomright',c('Wald','Logit'),col=cores[1:2],lty=1:2,pch=c(16,17),bty='n',cex=.7)
  plot(sub$n,sub$taxa_descarte,type='b',pch=16,col=cores[2],xlab='n',ylab='Discard rate',ylim=c(0,max(.01,sub$taxa_descarte)))
  dev.off()
  legendas <- c(legendas,paste0(tag,': true mu=',sub$mu_tau_verdadeiro[1],
    '; curves ordered by n=',paste(sub$n,collapse=', '),
    '; blue/orange/green/purple, line types 1/2/3/4. Black dashed reference: ',
    'true value (original scale and boxplots) or N(0,1) density (standardized scale).'))
}
write.csv(do.call(rbind,diagnosticos),file.path(saida,'Diagnostico_Densidades_QURS.csv'),row.names=FALSE)
writeLines(c('Results are conditional on replication acceptance.',
  paste('Standardization method:',centragem_padronizada),
  'SD is the empirical Monte Carlo standard deviation. The density comparison describes the shape of the standardized estimates.',
  legendas),file.path(saida,'Legendas_Figuras_QURS.txt'))
cat('Tables and figures saved to: ',saida,'\n')


