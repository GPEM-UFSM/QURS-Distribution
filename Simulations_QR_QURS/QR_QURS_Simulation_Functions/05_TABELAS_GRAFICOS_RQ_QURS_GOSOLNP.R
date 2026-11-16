# Generate LaTeX tables, performance plots, boxplots, and standardized densities.
# Run after script 03. Reporting requires a complete grid for at least one scenario.
# Pass a result directory name as the first Rscript argument.
args <- commandArgs(trailingOnly = TRUE)
nome_pasta_resultados <- if (length(args)) args[1L] else 'final_rq_qurs_gosolnp_v1_validas_R5000'
fazer_boxplots <- TRUE
fazer_densidades <- TRUE

.rqq_pasta<-local({
 a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
 fs<-Filter(function(x) !is.null(x),lapply(sys.frames(),function(x) x$ofile))
 f<-if(length(fs)) tail(fs,1)[[1]] else if(length(a)) sub('^--file=','',a[1]) else NULL
 if(is.null(f)) getwd() else dirname(normalizePath(f,mustWork=TRUE))
})
source(file.path(.rqq_pasta,'01_FUNCOES_RQ_QURS_GOSOLNP.R'),encoding='UTF-8')
source(file.path(.rqq_pasta,'02_MOTOR_MONTE_CARLO_RQ_QURS.R'),encoding='UTF-8')
raiz<-file.path(.rqq_pasta,'RESULTADOS_MONTE_CARLO',nome_pasta_resultados)
rqq_assert(file.exists(file.path(raiz,'Metadados_Estudo.rds')),'Study directory not found.')
meta<-readRDS(file.path(raiz,'Metadados_Estudo.rds'))
rqq_assert(identical(meta$versao,'RQ_QURS_gosolnp_v1'),
  'Pasta de outro modelo/versao: este codigo exige resultados RQ-QURS gosolnp v1.')
rqq_compilar_resultados(raiz)
controle<-read.csv(file.path(raiz,'Controle_Replicas_Geral.csv'),stringsAsFactors=FALSE)
resumo<-read.csv(file.path(raiz,'Resumo_Monte_Carlo_RQ_QURS_Geral.csv'),stringsAsFactors=FALSE)
ns<-sort(unique(meta$grade_completa$n));taus<-sort(unique(meta$grade_completa$tau));B<-meta$B
rqq_assert(length(taus)==5L && all(abs(taus-c(.1,.25,.5,.75,.9))<1e-12),
  'Este formato usa os cinco niveis tau do estudo.')
rqq_assert(length(ns)==4L,'Este formato de tabela usa quatro tamanhos; adapte multirow se mudar a grade.')
cenarios<-Filter(function(c) {
  esperado<-meta$grade_completa[meta$grade_completa$cenario==c,]
  linhas<-controle[controle$cenario==c,,drop=FALSE]
  ok<-nrow(linhas)==nrow(esperado) && !anyDuplicated(linhas$id_configuracao) &&
    setequal(linhas$id_configuracao,esperado$id_configuracao) &&
    all(linhas$completo) && all(linhas$replicas_validas==B)
  if(!ok) message('Skipping incomplete scenario: ',c)
  ok
},unique(meta$grade_completa$cenario))
rqq_assert(length(cenarios)>0,'No complete scenario is available for reporting.')
saida<-file.path(raiz,'TABELAS_E_GRAFICOS');dir.create(saida,showWarnings=FALSE)
rqq_fmt<-function(x,d=4L) {
 if(!is.finite(x)) return('--')
 if(abs(x)<.5*10^(-d)) x<-0
 paste0('$',sprintf(paste0('%.',d,'f'),x),'$')
}
rqq_nome_tex<-function(nm) paste0('$\\beta_{',sub('beta_','',nm,fixed=TRUE),'}$')
rqq_pdf<-function(arquivo,fun,width=9,height=7) {
 grDevices::pdf(arquivo,width=width,height=height)
 on.exit(grDevices::dev.off(),add=TRUE);fun()
}
todas<-character();medidas<-c( 'Mean'='media','Bias'='vies','SD'='desvio_padrao',
                             'MSE'='eqm','CP'='cobertura_wald')
for(c in cenarios) {
  sc<-meta$cenarios[[c]];tab<-resumo[resumo$cenario==c,,drop=FALSE]
  pasta<-file.path(saida,c);dir.create(pasta,showWarnings=FALSE)
  parametros<-unique(tab$parametro)
  for(param in parametros) {
    r<-tab[tab$parametro==param,,drop=FALSE]
    rqq_assert(nrow(r)==length(ns)*length(taus) && all(r$replicas==B),
      paste('Incomplete summary:',c,param))
    v<-unique(r$verdadeiro);rqq_assert(length(v)==1L,'Inconsistent true coefficient.')
    lines<-c('\\begin{table}[htbp]','\\centering',
      paste0('\\caption{Monte Carlo results for ',rqq_nome_tex(param),' in scenario ',c,
        ' of the QR-QURS model, with true value ',rqq_fmt(v),' e ',B,' valid replications per configuration.}'),
      paste0('\\label{tab:rqq_gosolnp_',tolower(c),'_',param,'}'),
      '\\resizebox{0.75\\textwidth}{!}{%','\\begin{tabular}{ccrrrrr}','\\toprule',
      'Measure & $n$ & $\\tau=0.10$ & $\\tau=0{,}25$ & $\\tau=0{,}50$ & $\\tau=0{,}75$ & $\\tau=0{,}90$ \\\\', '\\midrule')
    for(j in seq_along(medidas)) {
      if(j>1L) lines<-c(lines,'\\midrule')
      for(k in seq_along(ns)) {
        cells<-vapply(taus,function(t) {
          row<-r[abs(r$tau-t)<1e-12 & r$n==ns[k],,drop=FALSE]
          rqq_assert(nrow(row)==1L,'Missing or duplicate configuration.')
          rqq_fmt(row[[medidas[j]]])
        },'')
        lines<-c(lines,paste0(paste(c(if(k==1L) paste0('\\multirow{4}{*}{',names(medidas)[j],'}') else '',
          paste0('$',ns[k],'$'),cells),collapse=' & '),' \\\\'))
      }
    }
    lines<-c(lines,'\\bottomrule','\\end{tabular}','}',
      '\\begin{flushleft}\\scriptsize SD is the empirical standard deviation; MSE is the mean squared error; CP is the empirical coverage of 95\\% Wald confidence intervals on the coefficient scale. Results use accepted replications.\\end{flushleft}',
      '\\end{table}')
    writeLines(lines,file.path(pasta,paste0('Metricas_',param,'.tex')),useBytes=TRUE)
    todas<-c(todas,lines,'\\clearpage')
    # Plot bias, RMSE, and coverage by n and tau. Keep coverage on a common scale.
    rqq_pdf(file.path(pasta,paste0('Metricas_',param,'.pdf')),function() {
      par(mfrow=c(1,3),mar=c(4,4,3,1))
      for(metric in c('vies','raiz_eqm','cobertura_wald')) {
        ylim<-if(metric=='cobertura_wald') c(0,1) else range(c(0,r[[metric]]),finite=TRUE)
        if(diff(ylim)==0) ylim<-ylim+c(-.01,.01)
        plot(range(ns),ylim,type='n',xlab='n',ylab=c(vies='Bias', raiz_eqm='RMSE', cobertura_wald='Coverage')[metric],main=paste(c,param))
        for(t in seq_along(taus)) {
          rr<-r[abs(r$tau-taus[t])<1e-12,];rr<-rr[order(rr$n),]
          lines(rr$n,rr[[metric]],col=t+1L,lty=t,pch=14+t,type='b')
        }
        if(metric=='cobertura_wald') abline(h=.95,lty=3,col='gray40')
        if(metric=='vies') abline(h=0,lty=3,col='gray40')
        legend('topright',legend=paste('tau',taus),col=seq_along(taus)+1L,lty=seq_along(taus),pch=seq_along(taus)+14L,cex=.7,bty='n')
      }
    },width=13,height=4.5)
    for(tau in taus) {
      objetos<-lapply(ns,function(n) {
        f<-file.path(raiz,'CENARIOS',c,paste0('tau_',rqq_tag_tau(tau),'_n_',n),'Replicas_Validas.rds')
        o<-readRDS(f);rqq_assert(o$completo && length(o$replicas)==B,'Incomplete replication file.')
        vapply(o$replicas,function(x) x$coef[[param]],0.0)
      })
      tag<-paste0(param,'_tau_',rqq_tag_tau(tau))
      if(fazer_boxplots) rqq_pdf(file.path(pasta,paste0('Boxplot_',tag,'.pdf')),function() {
        boxplot(objetos,names=ns,xlab='n',ylab=param,main=paste(c,'tau=',tau),col='lightblue',outline=FALSE)
        abline(h=v,col='red3',lty=2)
      })
      if(fazer_densidades && B>=2L) rqq_pdf(file.path(pasta,paste0('Densidade_Padronizada_',tag,'.pdf')),function() {
        par(mfrow=c(2,2),mar=c(4,4,3,1));grid<-seq(-4,4,length.out=400)
        dens<-lapply(objetos,function(x) {
          if(!is.finite(sd(x)) || sd(x)==0) return(NULL)
          density((x-mean(x))/sd(x),from=-4,to=4,n=512)
        })
        ymax<-max(c(dnorm(0),unlist(lapply(dens,function(d) if(!is.null(d)) d$y))))*1.05
        for(k in seq_along(ns)) {
          plot(grid,dnorm(grid),type='l',lty=2,col='red3',ylim=c(0,ymax),
            xlab='(estimate - mean) / standard deviation',ylab='Density',main=paste(c,'tau=',tau,'n=',ns[k]))
          if(!is.null(dens[[k]])) lines(dens[[k]],col='navy',lwd=1.5)
          legend('topright',c('Monte Carlo','Normal(0,1)'),col=c('navy','red3'),lty=c(1,2),cex=.7,bty='n')
        }
      })
    }
  }
  rr<-controle[controle$cenario==c,];rr<-rr[order(rr$tau,rr$n),]
  lines<-c('\\begin{table}[htbp]','\\centering',
    paste0('\\caption{Replication counts for scenario ',c,' of the QR-QURS model (gosolnp).}'),
    paste0('\\label{tab:rqq_gosolnp_controle_',tolower(c),'}'),
    '\\begin{tabular}{ccrrrr}','\\toprule',
    '$\\tau$ & $n$ & Accepted & Attempts & Discarded & Discard rate (\\%) \\\\', '\\midrule')
  for(i in seq_len(nrow(rr))) lines<-c(lines,paste0(paste(c(rqq_fmt(rr$tau[i],2),rr$n[i],
    rr$replicas_validas[i],rr$tentativas[i],rr$replicas_descartadas[i],rqq_fmt(100*rr$taxa_descarte[i],2)),
    collapse=' & '),' \\\\'))
  lines<-c(lines,'\\bottomrule','\\end{tabular}',
    '\\begin{flushleft}\\scriptsize Discard rate equals 100 times the number discarded divided by the number attempted.\\end{flushleft}',
    '\\end{table}')
  writeLines(lines,file.path(pasta,'Controle_Replicas.tex'),useBytes=TRUE)
  todas<-c(todas,lines,'\\clearpage')
}
writeLines(todas,file.path(saida,'Todas_Tabelas_RQ_QURS.tex'),useBytes=TRUE)
writeLines(c('\\documentclass[11pt,a4paper]{article}','\\usepackage[utf8]{inputenc}',
  '\\usepackage[T1]{fontenc}','\\usepackage{booktabs,multirow,graphicx,amsmath}',
  '\\usepackage[margin=2cm]{geometry}','\\renewcommand{\\tablename}{Table}',
  '\\begin{document}',todas,'\\end{document}'),file.path(saida,'Documento_Tabelas_RQ_QURS.tex'),useBytes=TRUE)
message('Tables and plots saved to: ',saida)




