# use Lratio distance with Beta-binomial
DENDRO.simulation=function(kprob=NULL,lprob=NULL,filt=0,m=100,n=1000,
                           epi=0.001,RD=NULL,ref,k=NULL,subtype=1,rpt=100,
                           plot=TRUE){

  sampgenes=sample(seq(1,nrow(ref$X1)),n)
  sampcells=sample(seq(1,ncol(ref$X1)),m)
  sim=DENDROplan.simulation.helper2(kprob,lprob,epi,
                                    sampcells,sampgenes,
                                    RD,ref,k,subtype,verbose=TRUE)

  if(plot){
    plot(hclust(dist(t(sim$Z),'binary'),method='ward.D'),
         xlab='',ylab='',sub="",main='Example tree structure')
  }

  evals=replicate(rpt,DENDROplan.simulation.helper1(kprob=kprob,lprob=lprob,
                                                       filt,m,n,epi,RD,ref,k,
                                                       subtype,verbose=FALSE))

  res = matrix(NA,ncol=4,nrow=3)

  rownames(res)=c('CI_low','Mean','CI_up')
  colnames(res)=c('AdjustedRandIndex','Capture rate','Purity',
                   'Observation probability')

  res[1,c(1,2)]=pmax(apply(evals,1,function(x){
    sort(x[x>-999])[round(0.05*sum(x>-999))]
    }),0)[1:2]
  res[1,3]=pmax(apply(evals,1,function(x){
    sort(x[!is.na(x)])[round(0.05*sum(!is.na(x)))]
    }),0)[3]

  res[2,c(1,2)]=pmax(apply(evals,1,function(x){mean(x[x>-999])}),0)[1:2]
  res[2,3]=pmax(apply(evals,1,function(x){mean(x[!is.na(x)])}),0)[3]
  res[2,4]=apply(evals,1,function(x){sum(!is.infinite(x))/length(x)})[2]

  res[3,c(1,2)]=pmax(apply(evals,1,function(x){
    sort(x[x>-999])[round(0.95*sum(x>-999))]
  }),0)[1:2]
  res[3,3]=pmax(apply(evals,1,function(x){
    sort(x[!is.na(x)])[round(0.95*sum(!is.na(x)))]
  }),0)[3]

  res_tmp=res[,1:3]

  if(plot){
    g = ggplot2::ggplot(as.data.frame(t(rbind(res_tmp))),
                        ggplot2::aes(x = colnames(res_tmp), y = Mean)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_errorbar(ggplot2::aes(ymax = CI_up, ymin = CI_low)) +
      ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::ylim(0, 1)

    plot(g)
  }
  return(res)
}
