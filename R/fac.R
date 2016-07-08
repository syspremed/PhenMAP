fac <-
function(Y, X, minq=1, maxq=2, epsilon = 0.01, plot.BIC=FALSE, printout=TRUE) 
{
  V <- 5000  # maximum number of iterations
  N <- nrow(Y)
  p <- ncol(Y)
  L <- ncol(X)
  Covars <- cbind(X[,1],standardize(X[,-1,drop=FALSE]))
  Covars <- t(Covars)    # covariates
  
  #### storage of maxq estimates
  Sig_q <- matrix(NA, p,maxq)
  S2_q<-list()
  W_q <- list()
  U_q <- list()
  beta_q <- list()
  BIC <- rep(0, maxq)
  ll <- matrix(0, maxq, V)
  lla <- matrix(0, maxq, V)
  
  #### intializing the model and some prior values
  tcov<-t(Covars)
  icovtcov<-chol2inv(chol(Covars%*%tcov))
  tY<-t(Y)
  Ip<-diag(p)
  S <- (1/N)*(tY%*%Y)      # data covariance matrix
  temp <- eigen(S)         # eigen decomposition of the data covariance matrix 
  Vp <- 1                  # values for the prior covariance
  C2p <- 2                 # values for the prior covariance
  vo <- 1                  # values for the prior latent covariance
  Cp <-2                   # values for the prior latent covariance
  for (q in minq:maxq) {
    Iq<-diag(q)
    W<-temp$vec[,1:q]            # initializing the loadings matrix with eigen decomposition
    iSig<-Ip*(1/rep(abs((1/(p-q))*sum(temp$val[(q+1):p])),p)) # initializing the variance parameters with eigen decomposition
    u<-t(cmdscale(dist(Y),q))    # initializing the scores using MDS
    beta <- matrix(0, q, L)
    s2<-rep(0,q)
    for (i in 1:q){
      dat<-data.frame(u[i,],tcov[,-1,drop=FALSE])
      lm<-glm(dat, family = gaussian)    # initializing the regression coeffients with glm
      beta[i,]<-lm$coef 
      s2[i]<-var(lm$residuals)           # initializing the latent variance parameters with glm
    }
    iS2<-Iq*(1/s2)               
    del<-beta%*%Covars
    tol<-epsilon+1
    v<-0
    # EM algorithm
    while(tol>epsilon){
      v<-v+1
      ## LVs
      tws<-(t(W)%*%iSig)
      M_1<-chol2inv(chol( tws%*%W+iS2 ))
      u <- M_1%*%(tws%*%tY + iS2%*%del)
      tu<-t(u)
      ## Beta
      beta<-u%*%tcov%*%icovtcov
      del<-beta%*%Covars
      ## s2
      Sum_Euu <- (N*M_1) + (u%*%tu)
      MLES2 <- diag( Sum_Euu - 2*u%*%t(del) + del%*%t(del) )/N
      S2 <- (N*MLES2+Cp)/(N+vo+2)
      iS2<-Iq*(1/S2)
      ## loadings
      W <- (tY%*%tu)%*%chol2inv(chol(Sum_Euu))
      tW<-t(W)
      ## sig
      YEuW <- tY%*%tu%*%tW 
      MLESig <- diag( N*S - 2*YEuW + (W%*%Sum_Euu)%*%tW )/N
      Sig <- (N*MLESig+C2p)/(N+Vp+2)
      iSig<-Ip*(1/Sig)
      Den<-rep(NA,N)
      Sigma<-(W%*%((S2*Iq)%*%tW)+Sig*Ip)
      ## likelihood
      mumat<-W%*%del
      for(i in 1:N){
        Den[i]<-mvtnorm::dmvnorm(Y[i,],mumat[,i],Sigma,log=TRUE)
      }
      ll[q, v] <- sum(Den)
      plot(ll[q,1:v],col=2,xlab="EM iterations",ylab="log-likelihood",main=paste("model",q,sep="-")) 
      #if(v>1){if(ll[q,v-1]>ll[q,v]){stop("STOP: the likelihood is decreasing!")}}
      ## assess convergence
      converge <- Aitken(ll, lla, v, q, epsilon)
      tol <- converge[[1]]
      lla[q, v] <- converge[[2]]
      if (v == V) {
        if (printout == TRUE) {
          cat("Algorithm stopped for q = ", q, ". Maximum number of iterations exceeded.\n\n")
        }
        tol <- epsilon - 1
      }
    }# EM converged
    if (printout == TRUE) {
      cat("sFAC model:", q, "converged.\n\n")
    }
    params<-(p*q)-(0.5*q*(q-1))+p+q*(L+1)
    BIC[q]<-(2*ll[q,v]) - (params*log(N))
    U_q[[q]]<-u
    Sig_q[,q]<-Sig
    W_q[[q]]<-W
    beta_q[[q]] <- beta
    S2_q[[q]]<-S2
  }
  # Outputs of the optimal model via BIC
  qopt <- c(minq:maxq)[BIC[minq:maxq] == max(BIC[minq:maxq])]
  Uopt <- t(U_q[[qopt]])
  Wopt <- W_q[[qopt]]
  Sigopt <- Sig_q[,qopt]
  betaopt <- beta_q[[qopt]]
  S2opt <- S2_q[[qopt]]
  if (plot.BIC == TRUE) {
    plot(minq:maxq, BIC[minq:maxq], type = "b", xlab = "q", 
         ylab = "BIC", col.lab = "blue")
    abline(v = qopt, col = "red", lty = 2)
  }
  colnames(betaopt)<-colnames(X)
  rownames(betaopt)<-colnames(Wopt)<-colnames(Uopt)<-names(S2opt)<-paste("q",1:qopt,sep="_")
  names(Sigopt)<-rownames(Wopt)
  list(qopt = qopt, sig = Sigopt, factors = Uopt, loadings = Wopt, coefficients = betaopt, s2 = S2opt, covars=Covars, BIC = BIC[minq:maxq])
}
