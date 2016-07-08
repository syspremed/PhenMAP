mcmc.sampler <-
function(facfit, Y, X, w0, g0, mcmc, thinstep)
{
  covars<-t(X)
  L=nrow(covars)
  q=facfit$qopt
  n<-nrow(Y)
  p<-ncol(Y)
  S=mcmc
  nthin=1000
  thin=seq(1,S,by=thinstep)
  tY<-t(Y)
  In<-diag(n)
  Ip<-diag(p) 
  Iq<-diag(q)
  Il<-diag(L)
  rp<-rep(0,p)
  rq<-rep(0,q) 
  covtcov<-covars%*%t(covars)
  icovtcov<-chol2inv(chol( covtcov ))
  
  #### storage of S samples
  loadings.store<-array( NA, c(p,q,nthin) )
  scores.store<-array( NA, c(q,n,nthin) )
  sig.store<-matrix( NA, p, nthin )
  beta.store<-array( NA, c(q,L,nthin) )
  s2.store<-matrix( NA, q, nthin )
  
  #### initializing the model using sFAC, and setting prior values
  templateW<-W<-facfit$loadings   # initialize loading matrix using EM results
  a=b=w0   # priors for loadings matrix
  Sig<-facfit$sig; isig<-(1/Sig); iSig<-isig*Ip  # initialize variance parameters using EM results
  Sigo=nu0=.01; rS<-rep((nu0+n/2),p)  # priors for variance parameters
  Beta<-facfit$coeff  # initialize regression parameters using EM results
  g=g0; gtag<-g/(g+1)  # priors for regression parameters
  bcov<-gtag*icovtcov
  Hg<-(In - gtag*t(covars)%*%icovtcov%*%covars)
  s2<-facfit$s2; iS2<-(1/s2)*Iq  # initialize latent variance parameters using EM results
  s2o=nu0k=.01; rS2<-rep((nu0k+n/2),q)  # priors for latent variance parameters
    
  #### running the MCMC sampler
  set.seed(1)
  for(s in 1:S)
  {  
    if((round(s/500) - (s/500)) == 0){print(s)}         # Printing the number of iterations
    
    #### update scores u
    tWiS<-t(W)%*%iSig
    v.u<-chol2inv(chol(tWiS%*%W+iS2))
    m.u<-v.u%*%(tWiS%*%tY+iS2%*%Beta%*%covars)  
    U<-t(rmvnorm(n,rq,v.u))+m.u
    tU<-t(U)
    
    #### update loadings
    UU<-U%*%tU
    UY<-U%*%Y 
    for(j in 1:p)
    {
      w2<-(W[j,])^2
      theta<-rgamma(q,.5+a,.5*w2+b)
      v.l<-chol2inv( chol(UU*isig[j]+Iq*theta))
      m.l<-(v.l%*%UY[,j,drop=FALSE])*isig[j]
      W[j,]<-rmvnorm(1,m.l,v.l)  
    }#k
    
    #### update sig
    resid<-tY-W%*%U
    SSR<-diag(resid%*%t(resid)) 
    isig<-rgamma(p,rS,Sigo+.5*SSR) 
    iSig<-isig*Ip
    Sig<-1/isig
    
    #### update s2
    RSS<-diag(U%*%Hg%*%tU)
    is2<-rgamma(q,rS2,s2o+.5*RSS)
    iS2<-is2*Iq 
    s2<-1/is2
    
    #### update Beta
    ctU<-covars%*%tU
    for(k in 1:q)
    {
      v.b<-bcov*s2[k]
      m.b<-v.b%*%ctU[,k,drop=FALSE]*is2[k]
      Beta[k,]<-rmvnorm(1,m.b,v.b)       
    }#k   
    
    ## Storing parameter estimates and latent variables
    scores.store[,,s==thin] <- U
    loadings.store[,,s==thin] <- W 
    sig.store[,s==thin] <- Sig
    beta.store[,,s==thin] <- Beta
    s2.store[,s==thin] <- s2
  }#s  
  list(scores.store=scores.store, loadings.store=loadings.store, sig.store=sig.store, beta.store=beta.store, s2.store=s2.store, nthin=nthin)
}
