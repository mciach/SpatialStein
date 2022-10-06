###########load dependent libraries

library("Biobase")
library(Cardinal)
library(spam)
library(mvtnorm)
library("mclust")
#library(gridExtra)
library("ggplot2")
library(RColorBrewer)
################GMM
GMM<-function(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
{
  
  lod<-0
  x<-spectra(msset)[f,]
  ###leave out pixels with intensity lower than LOD and outliers
  if (length(x[x==min(x)])/ncol(msset)>0.01)
  {
    lod<-1
    x[x==min(x)]<-NA
    x2<-x[x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2])]
    x2<-x2[!is.na(x2)]
    
  } else
  {
    x2<-x[x<quantile(x)[3]+out*(quantile(x)[4]-quantile(x)[2])]
  }
  #####
  if (kprior==0)
  {
    gmm<-densityMclust(x2,modelNames="V")
    if (length(gmm$BIC)>2)
    {
      Di<-rep(NA,length(gmm$BIC)-1)
      for (i in 2:length(gmm$BIC))
      {
        Di[i-1]<-gmm$BIC[i]-gmm$BIC[i-1]
      }
      if (length(which(abs(Di)<3))!=0)
      {
        k<-min(which(abs(Di)<3))
      } else
      {
        k<-gmm$G
      }
    } else
    {
      k<-gmm$G
    }
    if (lod==0)
    {
      k<-min(k,kmax)
      kp<-k
    } else
    {
      k<-min(k,kmax-1)
      kp<-k+1
    }
  } else
  {
    k=kprior
    kp<-k
  }
  gmm<-densityMclust(x2,G=k,modelNames="V")
#  plot(gmm,what="density",x2,breaks=50)
  aa<-rep(0,ncol(msset))
  
  aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-apply(gmm$z,1, function (x) which(x==max(x)))
  msset$gmm<-aa
#  image(msset, formula = gmm~x*y,main=paste0("feature",f))
  k<-kp
  return(list(gmm,k))
}


#########################DGMM with different radii
K_DGMM<-function(msset=msset,gmm=gmm,f=f,k=k,step=1e7,initialization="km",r_max=3)
{
  kr<-rep(0,r_max)
  for (radius in 1:r_max)
  {
    ##################neighboring matrix
    w<-list()
    coords<-coord(msset)
    for (i in 1:ncol(msset))
    {
      w[[i]] <- 
        (abs(as.numeric(coords[i,'y']) - as.numeric(coords[,'y'])) <= radius) & (abs(as.numeric(coords[i,"x"]) - as.numeric(coords[,"x"])) <= radius)
      w[[i]]<-which(w[[i]]==T)
      w[[i]]<-w[[i]][w[[i]]!=i]
    }
    
    
    rmlist<-which(is.na(w)==T)
    if (length(rmlist)!=0)
    {
      w<-w[[-rmlist]]
    }
    
    ###################################fit DGMM candidate
    ##########ion intensities
    int<-spectra(msset)[f,]
    if (length(rmlist)!=0)
    {
      int<-int[-rmlist]
    }
    x<-int
    #######number of pixels
    N=length(x)
    
    K<-k
    g<-k
    ###############initialize using k-means
    
    if (initialization=="km")
    {
      km<-kmeans(int,centers =k)
      mu<-km$centers
      sigma<-(mu*0.2)^2
 
      
    }
    
    ##############initialize using Gaussian Mixture Model
    
    if (initialization=="gmm")
    {
      if (gmm$G != k)
      {
        mu[1]<-0.1
        sigma[1]<-0.004
        mu[2:k]<-gmm$parameters$mean
        sigma[2:k]<-gmm$parameters$variance$sigmasq
      } else
      {
        mu<-gmm$parameters$mean
        sigma<-gmm$parameters$variance$sigmasq

      }
      
    }
    
    ############initialize alpha in Dirichlet process
    alpha=rep(1,g);
    ############initialize beta in PI
    beta=1;
    
    ###########step size
    eta<-min(mu)/step
    ##########differentials of mu, sigma and alpha
    dmu<-rep(1,g)
    dsg<-rep(1,g)
    dalpha<-rep(1,g)
    #########posterior probability
    y<-matrix(0, nrow=N, ncol=K)
    #########prior probability
    PI<-matrix(1/K, nrow=N, ncol=K)
    logPI<-matrix(1/K, nrow=N, ncol=K)
    #########P(x|mu, sigma)
    px<-matrix(0, nrow=N, ncol=K)
    logpx<-matrix(0, nrow=N, ncol=K)
    iteration=100
    #########trace
    mutrace<-matrix(0,ncol=K,nrow=iteration)
    sigtrace<-matrix(0,ncol=K,nrow=iteration)
    alphatrace<-matrix(0,ncol=K,nrow=iteration)
    betatrace<-matrix(0,ncol=1,nrow=iteration)
    #########negative loglikelihood
    loglik<-rep(0,iteration)
    #########initialize P(x|mu,sigma)
    for (j in 1:K)
    {
      px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
    }
    
    for (j in 1:K)
    {
      logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
    }
    ######### initialize posterior probability
    y<-px*PI/rowSums(px*PI)
    y[is.na(y)==TRUE]<-1/k
    y[y==0]<-1e-200
    ybar<-y
    for (i in 1:iteration)
    {
      
      
      ############average posterior probability
      for (j in 1:ncol(msset))
      {
        ybar[j,]<-colMeans(y[w[[j]],])
      }
      
      ybar[ybar==0]<-1e-100
      
      #############negative loglikelihod
      loglik[i]<--sum(log(rowSums(t(t((ybar)^beta)*alpha^2)/rowSums(t(t((ybar)^beta)*alpha^2))*px)))
      logybar<-log(ybar)
      for ( j in 1:K)
      {
        logPI[,j]<-2*log(abs(alpha[j]))+beta*logybar[,j]
      }
      logPI<-logPI-rowMin(logPI)
      
      for ( j in 1:K)
      {
        PI[,j]<-alpha[j]^2*ybar[,j]^beta
      }
      PI[PI==Inf]<-1e100
      PI<-PI/rowSums(PI)
      
      PI[PI==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
      }
      ##logp(x|mu, sigma)
      for (j in 1:K)
      {
        logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
      }
      ##posterior
      
      
      y<-px*PI/rowSums(px*PI)
      y[is.na(y)==TRUE]<-1/k
      y[y==0]<-1e-100
      for ( j in 1:K)
      {
        dmu[j]<-sum(y[,j]*1/sigma[j]*(mu[j]-x))
        dsg[j]<-1/2*sum(y[,j]*(1/sigma[j]-(x-mu[j])^2/(sigma[j]^2)))
        dalpha_p<-y*(ybar[,j])^beta/rowSums(t(t((ybar)^beta)*alpha^2))
        dalpha_p[is.na(dalpha_p)]<-1
        dalpha[j]<--sum(2*y[,j]/alpha[j])+2*alpha[j]*sum(dalpha_p)
        
      }
      
      dbeta=sum(y*(-log(ybar)+rowSums(t(t((ybar)^beta)*alpha^2)*log(ybar))/rowSums(t(t((ybar)^beta)*alpha^2))))
      mu<-mu-eta*dmu
      sigma<-sigma-eta*dsg
      sigma[sigma<=0]<-0.006327605
      alpha<-alpha-eta*dalpha
      beta<-beta-eta*dbeta
      beta<-max(beta,0)
      beta<-min(beta,10)
      mutrace[i,]<-mu
      sigtrace[i,]<-sigma
      alphatrace[i,]<-alpha
      betatrace[i,]<-beta
      
    }
    xx<-rep(1,ncol(msset))
    if (length(rmlist)!=0)
    {
      xx[-rmlist]<-apply(ybar,1, function (x) which(x==max(x))) 
    } else
    {
      xx<-apply(ybar,1, function (x) which(x==max(x)))
    }
    
    msset$dgmm<-xx
#    image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F)
    kprim<-length(unique(msset$dgmm))
    L<-unique(msset$dgmm)
    for (i in L)
    {
      if (length(msset$dgmm[msset$dgmm==i])/ncol(msset)<0.01)
      {
        kprim<-kprim-1
        msset$dgmm[msset$dgmm==i]<-NA
      }
      
    }
    print(kprim)
    print(min(loglik))
    print(mu)
    seg<-unique(msset$dgmm[!is.na(msset$dgmm)])
    seg2<-seg
    for (i in 1:length(seg[!is.na(seg)]))
    {
      j<-i+1;
      while(j<length(seg[!is.na(seg)]))
      {
        if ((mu[seg[i]]-mu[seg[j]])/sqrt(sigma[seg[i]]/2+sigma[seg[j]]/2)<1)
        {
          seg2[j]<-NA
        }
        j=j+1
      }
    }
    kr[radius]<-length(seg2[!is.na(seg)])
  }
  kr<-Mode(kr)
  print(kr)
  return(kr)
}

############################spatial DGMM fit

DGMM<-function(msset=msset,w,w2=w2,k,f,sp_ratio=4,step=1e5, iteration=1000,Annl=0, sig=0.2,initialization="km")
{
  
  ############# remove isolated pixels
  
  rmlist<-NULL
  
  
  if (length(rmlist)!=0)
  {
    w<-w[-rmlist,-rmlist]
  }
  
  ###################################fit DGMM
  int<-spectra(msset)[f,]
  if (length(rmlist)!=0)
  {
    int<-int[-rmlist]
  }
  #Annl=0
  tt<-1
  
  x<-int
  N=length(x)
  ###############initialize using k-means
  
  if (initialization=="km")
  {
    km<-kmeans(int,centers =k)
    mu<-sort(km$centers, decreasing = F)
    sigma<-(mu*sig)^2
    sigma[sigma<0.0006]<-0.0006
    print(mu)
    print(sigma)
  }
  
  ##############initialize using Gaussian Mixture Model
  
  if (initialization=="gmm")
  {
    lod<-0
    x<-int
    ###leave out pixels with intensity lower than LOD and outliers
    if (length(x[x==min(x)])/ncol(msset)>0.01)
    {
      lod<-1
      x[x==min(x)]<-NA
      x2<-x[x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2])]
      x2<-x2[!is.na(x2)]
      
    } else
    {
      x2<-x[x<quantile(x)[3]+out*(quantile(x)[4]-quantile(x)[2])]
    }
    gmm<-densityMclust(x2,G=k,modelNames="V")

    mu<-gmm$parameters$mean
    sigma<-gmm$parameters$variance$sigmasq
    print(mu)
    print(sigma)
    #  plot(gmm,what="density",x2,breaks=50)
    aa<-rep(0,ncol(msset))
    
    aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-apply(gmm$z,1, function (x) which(x==max(x)))
    msset$gmm<-aa
    print(image(msset, formula = gmm~x*y,main=paste0("feature",f)))
  }
  
  ############initialize alpha in Dirichlet process
  
  alpha<-rep(1,k)
  ############initialize beta in PI
  beta=5;
  K<-k
  #########step size
  eta<-min(mu)/step
  ###########differentials of mu, sigma and alpha
  dmu<-rep(1,k)
  dsg<-rep(1,k)
  dalpha<-rep(1,k)
  ###########posterior probability
  y<-matrix(0, nrow=N, ncol=K)
  #########prior probability
  PI<-matrix(1/K, nrow=N, ncol=K)
  logPI<-matrix(1/K, nrow=N, ncol=K)
  #########P(x|mu, sigma)
  px<-matrix(0, nrow=N, ncol=K)
  logpx<-matrix(0, nrow=N, ncol=K)
  #iteration=100
  #########trace
  mutrace<-matrix(0,ncol=k,nrow=iteration)
  sigtrace<-matrix(0,ncol=k,nrow=iteration)
  alphatrace<-matrix(0,ncol=K,nrow=iteration)
  betatrace<-matrix(0,ncol=1,nrow=iteration)
  #########negative loglikelihood
  loglik<-rep(0,iteration)
  PI_p<-matrix(0, nrow=N, ncol=K)
  px_p<-matrix(0, nrow=N, ncol=K)
  #########initialize P(x|mu,sigma)
  for (j in 1:K)
  {
    px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
  }
  ##logp(x|\theta)
  for (j in 1:K)
  {
    logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
  }
  ######## initialize posterior probability
  
  y<-px*PI/rowSums(px*PI)
  y[is.na(y)==TRUE]<-runif(length(y[is.na(y)==TRUE]),0,1)
  ##########handling data out of storage range
  y[y==0]<-1e-100
  
  ybar<-y
  for (i in 1:iteration)
  {
    ############average posterior probability
    for (j in 1:ncol(msset))
    {
      ybar[j,]<-w2[[j]]%*%y[w[[j]],]/sum(w2[[j]]%*%y[w[[j]],])
    }
    
    ybar[ybar==0]<-1e-100
    ybar[is.na(ybar)==TRUE]<-runif(length(ybar[is.na(ybar)==TRUE]),0,1)
    #############negative loglikelihod
    loglik[i]<--sum(log(rowSums(t(t((ybar)^beta)*alpha^2)/rowSums(t(t((ybar)^beta)*alpha^2))*px)))
    logybar<-log(ybar)
    
    
    for ( j in 1:K)
    {
      PI[,j]<-alpha[j]^2*ybar[,j]^beta
    }
    PI[PI==Inf]<-1e100
    PI<-PI/rowSums(PI)
    
    PI[PI==0]<-1e-100 
    ##p(x|mu, sigma)
    for (j in 1:K)
    {
      px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
    }
    ##logp(x|mu, sigma)1
    for (j in 1:K)
    {
      logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
    }
    ##posterior
    px<-px/rowSums(px)
    
    y<-px*PI/rowSums(px*PI)
    y[is.na(y)==TRUE]<-runif(length(y[is.na(y)==TRUE]),0,1)
    y[y==0]<-1e-100
    
    ###################calculate differential
    for ( j in 1:K)
    {
      dmu[j]<-sum(y[,j]*1/sigma[j]*(mu[j]-x))
      dsg[j]<-1/2*sum(y[,j]*(1/sigma[j]-(x-mu[j])^2/(sigma[j]^2)))
      dalpha_p<-y*(ybar[,j])^beta/rowSums(t(t((ybar)^beta)*alpha^2))
      dalpha_p[is.na(dalpha_p)]<-1
      dalpha[j]<--sum(2*y[,j]/alpha[j])+2*alpha[j]*sum(dalpha_p)
      
    }
    
    dbeta=sum(y*(-log(ybar)+rowSums(t(t((ybar)^beta)*alpha^2)*log(ybar))/rowSums(t(t((ybar)^beta)*alpha^2))))
    ##################updata parameters
    mu<-mu-eta*dmu
    for (mu_ind in 1:k)
    {
      mu[mu_ind]<-max(mu[mu_ind],0)
    }
    sigma<-sigma-eta*dsg
    sigma[sigma<=0]<-0.006327605
    alpha<-alpha-eta*dalpha
    beta<-beta-eta*dbeta
    beta<-max(beta,0)
    beta<-min(beta,10)
    
    ############################simulated annealing
    if (Annl==TRUE)
    {
      ####randomly change first element of mu
      mu_p<-c(runif(1, min=max(mu[1]-0.2*mu[1]*tt,0), max=mu[1]+0.2*mu[1]*tt), mu[2:k])
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2
      }
      PI_p<-PI/rowSums(PI)
      
      PI_p[PI_p==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu_p[j])^2/2/sigma[j])
      }
      y_p<-PI_p*px_p/rowSums(PI_p*px_p)
      y_p[y_p==0]<-1e-100
      ybar_p<-y_p
      for (j in 1:ncol(msset))
      {
        ybar[j,]<-w2[[j]]%*%y[w[[j]],]/sum(w2[[j]]%*%y[w[[j]],])
      }
      ybar_p[ybar_p==0]<-1e-100
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2*ybar_p[,j]^beta
      }
      PI_p<-PI_p/rowSums(PI_p)
      
      PI_p[PI_p==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu_p[j])^2/2/sigma[j])
      }
      
      loglik1<--sum(log(rowSums(t(t((ybar_p)^beta)*alpha^2)/rowSums(t(t((ybar_p)^beta)*alpha^2))*px_p)))
      #######################
      for (j in 1:ncol(msset))
      {
        ybar_p[j,]<-w2[[j]]%*%y[w[[j]],]/sum(w2[[j]]%*%y[w[[j]],])
      }
      ybar_p[ybar_p==0]<-1e-100
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2*ybar[,j]^beta
      }
      PI_p<-PI/rowSums(PI)
      
      PI_p[PI_p==0]<-1e-100 
      
      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
      }
      
      loglik2<--sum(log(rowSums(t(t((ybar_p)^beta)*alpha^2)/rowSums(t(t((ybar_p)^beta)*alpha^2))*px_p)))
      
      
      if (loglik1<loglik2)
      {
        mu<-mu_p
        y<-y_p
      }
    }
    ###################################################################
    
    mutrace[i,]<-mu
    sigtrace[i,]<-sigma
    alphatrace[i,]<-alpha
    betatrace[i,]<-beta
    #########cooling
    tt<-tt-1/iteration
    
    
  }
  
  
  
  xx<-rep(1,ncol(msset))
  if (length(rmlist)!=0)
  {
    xx[-rmlist]<-apply(y,1, function (x) which(x==max(x))) 
  } else
  {
    xx<-apply(y,1, function (x) which(x==max(x)))
  }
  
  msset$dgmm<-xx
 print(image(msset, formula = dgmm~x*y,asp=sp_ratio, main=paste0(f,"y")))
  
  xx<-rep(1,ncol(msset))
  if (length(rmlist)!=0)
  {
    xx[-rmlist]<-apply(ybar,1, function (x) which(x==max(x))) 
  } else
  {
    xx<-apply(y,1, function (x) which(x==max(x)))
  }
  
  msset$dgmm<-xx
#  image(msset, formula = dgmm~x*y,asp=sp_ratio, main=paste0(f,"ybar"))
  
  return(list(mu,sigma,alpha,beta,mutrace,sigtrace,alphatrace,betatrace,loglik,y, xx))
}


W_matrix<-function(msset,sp_ratio=2,radius=1,sigma=1e4)
{
  ##################weight matrix
  w<-list()
  radius<-1 
  coords<-coord(msset)
  for (i in 1:ncol(msset))
  {
    w[[i]] <- 
      (abs(as.numeric(coords[i,'y']) - as.numeric(coords[,'y'])) <= radius) & (abs(as.numeric(coords[i,"x"]) - as.numeric(coords[,"x"])) <= radius)
    w[[i]]<-which(w[[i]]==T)
    w[[i]]<-w[[i]][w[[i]]!=i]
  }
  
  
  #################similarity score
  w2<-w
  rccnorep<-msset
  for (i in 1:ncol(msset))
  {
    for (k in 1:length(w[[i]]))
    {
      j<-w[[i]][k]
      w2[[i]][k]<-exp(-((coord(rccnorep)$x[i]-coord(rccnorep)$x[j])^2+sp_ratio*(coord(rccnorep)$y[i]-coord(rccnorep)$y[j])^2))*exp(-t(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])%*%(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])/sigma)
      
    }
  }
  return(list(w,w2))
}


Mode<-function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



GMM<-function(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
{
  
  lod<-0
  x<-spectra(msset)[f,]
  #x = x - min(x)
  #x = log(x+1)
  ###leave out pixels with intensity lower than LOD and outliers
  if (length(x[x==min(x)])/ncol(msset)>0.01)
  {
    lod<-1
    x[x==min(x)]<-NA
    x2<-x[x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2])]
    x2<-x2[!is.na(x2)]
    
  } else
  {
    x2<-x[x<quantile(x)[3]+out*(quantile(x)[4]-quantile(x)[2])]
  }
  #####
  if (kprior==0)
  {
    gmm<-densityMclust(x2,modelNames="V")
    id = 1:length(gmm$BIC)
    id = id[is.na(gmm$BIC)==F]
    if (length(gmm$BIC[is.na(gmm$BIC)==F])>2)
    {
      bic = gmm$BIC[is.na(gmm$BIC)==F]
      ncomp <- length(bic)
      b<-c()
      b[1]<-ncomp-1
      b[2]<-bic[ncomp]-bic[1]
      max_d=0
      max_c=1
      for (i in 1:ncomp)
      {
        
        p<-c()
        p[1] = i-1
        p[2] = bic[i]-bic[1]
        pe<-sum(p*b)/sqrt(sum(b^2))*b/sqrt(sum(b^2))
        d<-sum((p-pe)^2)
        print(d)
        if (d > max_d)
        {
          
          max_d = d
          max_c = id[i]
        }
      }
    }else
    {
      max_c = which(gmm$BIC == max(gmm$BIC,na.rm = T))
    }
  }else
  {
    max_c = kprior
  }
  gmm<-densityMclust(x2,G=max_c,modelNames="V")
  #  plot(gmm,what="density",x2,breaks=50)
  aa<-rep(0,ncol(msset))
  
  aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-apply(gmm$z,1, function (x) which(x==max(x)))
  msset$gmm<-aa
  aa[aa==0]<-1
  #print(image(msset, formula = gmm~x*y,main=paste0("feature",f)))
  return(list(gmm,max_c,aa))
}

