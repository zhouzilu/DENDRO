
DENDROplan.simulation.helper1=function(kprob,lprob,
                           filt,m,n,
                           epi,RD,ref,k,subtype,verbose){
  sampgenes=sample(seq(1,nrow(ref$X1)),n)
  sampcells=sample(seq(1,ncol(ref$X1)),m)
  sim=DENDROplan.simulation.helper2(kprob,lprob,epi,
                                   sampcells,sampgenes,
                                   RD,ref,k,subtype,verbose)
  return(DENDROplan.simulation.helper3(sim$X,sim$N,sim$Z,sim$clades,filt))
}


# Simulation of N, X, Z matrix in a 3 clades tree. On top of SNV_simulation 4,
# we could handle different read depth with sample than stupid multiplication
# We also apply the newest distance calculation method (incoporate
# beta-binomial model)
DENDROplan.simulation.helper2=function(kprob=NULL,lprob=NULL,
                                      epi=0.001,sampcells,
                                      sampgenes,RD,ref,k, subtype,verbose){
  #*8464/200 #based on http://onlinelibrary.wiley.com/doi/10.1111/j.1755-0998.2011.03024.x/full
  X1.samp=ref$X1[sampgenes,sampcells]
  X2.samp=ref$X2[sampgenes,sampcells]
  # print(mean(CAST_rc.samp,na.rm=T))
  # genenames.samp=genenames[sampgenes]
  # cellnames.samp=cellnames[sampcells]
  # cell.stage.samp=cell.stage[sampcells]
  if(!is.null(RD)){
    samp_RD=mean(ref$X1+ref$X2,na.rm=TRUE)
    RD_ratio=RD/samp_RD
    X1.samp=matrix(rbinom(length(as.vector(X1.samp)),
                          as.vector(X1.samp),RD_ratio),nrow=nrow(X1.samp))
    X2.samp=matrix(rbinom(length(as.vector(X2.samp)),
                          as.vector(X2.samp),RD_ratio),nrow=nrow(X2.samp))
  }

  # k is the number of clade. Default 3.
  if(!is.null(kprob)){
    k=length(kprob)
  }else if (is.null(k)){
    k=3
  }
  if(verbose){
    cat('There are total ',k,' clades in our simulation \n')
  }


  m=length(sampcells)
  n=length(sampgenes)
  # cat('in SNV_simulation4, m:',m,' n:',n,'\n')
  if(is.null(kprob)){
    clades=sample(seq(1,k),length(sampcells),replace=T)
  }
  else{
    clades=sample(seq(1,k),length(sampcells),prob=kprob,replace=T)
  }

  if(is.null(lprob)){
    l=2*k-2
    cladegene=sample(seq(1,l),length(sampgenes),replace=T)
  }
  else{
    l=length(lprob)
    if (l<2*k-2){
      cat('Error with lprob, lprob should have length = 2*length of kprob - 2')
    }
    cladegene=sample(seq(1,l),length(sampgenes), prob=lprob,replace=T)
  }

  if(k==2){
    N=X=Z=matrix(NA,ncol=m,nrow=n)
    for(i in 1:n){
      if (sample(c('X1','X2'),1)=='X1'){
        Ncga=X1.samp[i,]
        Ncgb=X2.samp[i,]
      }
      else{
        Ncga=X2.samp[i,]
        Ncgb=X1.samp[i,]
      }
      N[i,]=Ncga+Ncgb
      if(cladegene[i]==1){
        for(j in 1:m){
          if (clades[j]==1){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
      else if (cladegene[i]==2){
        for(j in 1:m){
          if (clades[j]==2){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }else{
        for(j in 1:m){
          X[i,j]=rbinom(1,N[i,j],epi)
          Z[i,j]=0
        }
      }
    }
  }

  if(k==3){
    N=X=Z=matrix(NA,ncol=m,nrow=n)
    for(i in 1:n){
      if (sample(c('X1','X2'),1)=='X1'){
        Ncga=X1.samp[i,]
        Ncgb=X2.samp[i,]
      }
      else{
        Ncga=X2.samp[i,]
        Ncgb=X1.samp[i,]
      }
      N[i,]=Ncga+Ncgb
      if(cladegene[i]==1){
        for(j in 1:m){
          if (clades[j]==1){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
      else if (cladegene[i]==2){
        for(j in 1:m){
          if (clades[j]==2){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
      else if (cladegene[i]==3){
        for(j in 1:m){
          if (clades[j]==3){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
      else if (cladegene[i]==4){
        for(j in 1:m){
          if (clades[j]==2|clades[j]==3){
            X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
            Z[i,j]=1
          }
          else{
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
      else{
        for(j in 1:m){
          X[i,j]=rbinom(1,N[i,j],epi)
          Z[i,j]=0
        }
      }
    }
  }

  if(k==4){
    if(subtype==1){
      N=X=Z=matrix(NA,ncol=m,nrow=n)
      for(i in 1:n){
        if (sample(c('X1','X2'),1)=='X1'){
          Ncga=X1.samp[i,]
          Ncgb=X2.samp[i,]
        }
        else{
          Ncga=X2.samp[i,]
          Ncgb=X1.samp[i,]
        }
        N[i,]=Ncga+Ncgb
        if(cladegene[i]==1){
          for(j in 1:m){
            if (clades[j]==1){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==2){
          for(j in 1:m){
            if (clades[j]==2){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==3){
          for(j in 1:m){
            if (clades[j]==3){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==4){
          for(j in 1:m){
            if (clades[j]==4){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==5){
          for(j in 1:m){
            if (clades[j]==1|clades[j]==2){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==6){
          for(j in 1:m){
            if (clades[j]==3|clades[j]==4){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else{
          for(j in 1:m){
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
    }else{
      N=X=Z=matrix(NA,ncol=m,nrow=n)
      for(i in 1:n){
        if (sample(c('X1','X2'),1)=='X1'){
          Ncga=X1.samp[i,]
          Ncgb=X2.samp[i,]
        }
        else{
          Ncga=X2.samp[i,]
          Ncgb=X1.samp[i,]
        }
        N[i,]=Ncga+Ncgb
        if(cladegene[i]==1){
          for(j in 1:m){
            if (clades[j]==1){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==2){
          for(j in 1:m){
            if (clades[j]==2){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==3){
          for(j in 1:m){
            if (clades[j]==3){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==4){
          for(j in 1:m){
            if (clades[j]==4){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==5){
          for(j in 1:m){
            if (clades[j]==3|clades[j]==4){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else if (cladegene[i]==6){
          for(j in 1:m){
            if (clades[j]==2|clades[j]==3|clades[j]==4){
              X[i,j]=rbinom(1,Ncga[j],1-epi)+rbinom(1,Ncgb[j],epi)
              Z[i,j]=1
            }
            else{
              X[i,j]=rbinom(1,N[i,j],epi)
              Z[i,j]=0
            }
          }
        }
        else{
          for(j in 1:m){
            X[i,j]=rbinom(1,N[i,j],epi)
            Z[i,j]=0
          }
        }
      }
    }
  }
  return(list(N=N,X=X,Z=Z,clades=clades))
}



DENDROplan.simulation.helper3=function(X,N,Z,clades,filt,epi=0.001){
  stat_lr = tryCatch({
    # print(epi)
    # filter out genes that there are no mutations
    memb_true=clades
    num_memb_true=cumsum(table(clades))
    k=length(unique(clades))

    r_filt=rowSums(X>0,na.rm = T)>filt
    X_filt=X[r_filt,]
    N_filt=N[r_filt,]
    Z_filt=Z[r_filt,]
    ord=order(memb_true)
    N_filt_sort=N_filt[,ord]
    X_filt_sort=X_filt[,ord]
    Z_filt_sort=Z_filt[,ord]
    memb_true=memb_true[ord]
    d<-DENDRO.dist(X_filt_sort,N_filt_sort,Z_filt_sort,epi,show.progress=FALSE)
    # d=d/sd(d,na.rm=T)
    # d[which(is.na(d),arr.ind=TRUE)]=max(d,na.rm=TRUE)+1
    hc=hclust(d,method='ward.D')
    memb_pred=cutree(hc, k = k)
    stat1=mclust::adjustedRandIndex(memb_true, memb_pred)
    stat2=capture_rate(memb_true,memb_pred,k)
    stat3=purity(memb_true,memb_pred,k)
    c(stat1,stat2,stat3)
  },error=function(e){
    c(NA,NA,NA)
  })
  return(stat_lr)
}

capture_rate=function(memb_true,memb_pred,k){
  return(max(table(memb_pred[memb_true==k]))/sum(memb_true==k))
}

purity=function(memb_true,memb_pred,k){
  g=as.numeric(names(sort(table(memb_pred[memb_true==k]),decreasing=TRUE))[1])
  return(sum(memb_pred==g & memb_true==k)/sum(memb_pred==g))
}

