sBFAC <-
function(Y,X,fmin=1,fmax=2,gmin=1,gmax=2,w0=1,g0=nrow(Y),mcmc=20000)
{
  if (missing(Y)) {
    stop("Expression data are required to fit the model.\n")
  }
  if (missing(X)) {
    stop("Covariates for this model.\n ")
  }
  if (nrow(Y) != nrow(X)) {
    stop("Expression data and covariate data should have the same number of rows.\n")
  }
  if (missing(fmin)) {
    fmin <- 1
  }
  if (missing(fmax)) {
    fmax <- 2
  }
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(fmin)) {
    stop("fmin should be a whole number.\n")
  }
  if (!is.wholenumber(fmax)) {
    stop("fmax should be a whole number.\n")
  }
  if (fmin > fmax) {
    stop("fmin can not be greater than fmax\n")
  }
  if (fmax > ncol(Y)) {
    stop("fmax can not be greater than the number of variables.\n")
  }
  if (fmax > 10) {
    cat("Warning! Model fitting may become very slow for fmax > 10.\n\n")
  }
  if (missing(mcmc)) {
    mcmc <- 20000
  }
  if (mcmc < 1000) {
    stop("mcmc can not be less than 1000.\n")
  }
  if (mcmc < 1000) {
    stop("mcmc can not be less than 1000.\n")
  }
  if ((mcmc%%1000)!=0) {
    stop("mcmc should be devisible by 1000.\n")
  }
  if (mcmc < 10000) {
    cat("Warning! Poor convergence for mcmc less than 10000 .\n\n")
  }
  if (gmin<1) {
    stop("gmin can not be less than 1.\n")
  }
  if (!is.wholenumber(gmin)) {
    stop("gmin should be a whole number.\n")
  }
  if (gmax>.2*nrow(Y)) {
    stop("gmax can not more than a fifth of the samples.\n")
  }
  if (!is.wholenumber(gmax)) {
    stop("gmax should be a whole number.\n")
  }
  if ((w0<.01)|(w0>10)) {
    stop("w0 should not be more than 10 and less than 0.01.\n")
  }
  if (gmax<gmin) {
    stop("gmax can not less than gmin.\n")
  }
  if (g0<1) {
    stop("g0 can not be less than 1.\n")
  }
  ####################################
  ### installing missing R packages 
  installed <- installed.packages()[,1]                              
  required <- c("mvtnorm", "mclust", "MetabolAnalyze", "stats")
  toinstall <- required[!(required %in% installed)]
  if(length(toinstall) != 0){
    install.packages(toinstall)
   }
  lapply(required, require, character.only = TRUE)
  ###

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  p<-ncol(Y)  # number of features
  n<-nrow(Y)  # number of samples
  L<-ncol(X)  # number of covariates
  thinstep=mcmc/1000  # thinning steps
  burn=200    # burn in period 
  ################################################################################################
  #### (1) Select the optimal number of latent variables by Fitting several sBFAC via EM algorithm
  ####     varying the number of latent variables from fmin to fmax.
  par(mfrow=c(3,3))
  Y <- sweep(Y,2,colMeans(Y),"-")         # centering data                         
  res1<-fac(Y,X,fmin,fmax,plot.BIC=T)     
  plot(res1[[3]][,1:2],col=1,xlab="factor:1",ylab="factor:2")
  
  #########################
  #### (11) sBFAC MCMC algorithm
  res2<-mcmc.sampler(res1,Y,X,w0,g0,mcmc,thinstep)
  print("finished running MCMC sampler!!")
  
  ########################################################################################
  #### Model identification via Procrustes Post Processing: Rotate W to match W at burn-in
  #### then use the ration matrix on U and Beta
  q<-res1$qopt
  Iq<-diag(q)
  nthin<-res2$nthin
  u.store<-res2$scores.store
  w.store<-res2$loadings.store
  b.store<-res2$beta.store
  rot.loadings.store <- array(NA, c(p,q,nthin))
  rot.scores.store <- array(NA, c(q,n,nthin))
  rot.beta.store <- array(NA, c(q,L,nthin))
  rot.s2.store <- matrix(NA, q,nthin)
  for(i in 1:nthin)
  {
    if(q==1){
      ltemp<-matrix(w.store[,,i])  
      target<-matrix(res1$loadings)
    }else{
      ltemp<-w.store[,,i]     
      target<-res1$loadings
    }#q=1
    res <- MCMCpack:::procrustes(ltemp, target, translation=FALSE, dilation=FALSE)
    rot.loadings.store[,,i] <- res$X.new
    rot.scores.store[,,i] <- t(res$R)%*%u.store[,,i]
    rot.beta.store[,,i] <- t(res$R)%*%b.store[,,i]
    rot.s2.store[,i] <- res2$s2.store[,i]
  }#i
  
  ####################################
  #### Diagnostics plots: trace plots
  par(mfrow=c(3,3), mar=c(4,3.7,1.5,0.5), oma=c(0,0,0,0),mgp=c(2.5,1,0))
  matplot(t(rot.scores.store[sample(1:q,1),sample(1:n,ceiling(.25*n)),burn:nthin]), type="l",main="Factors trace plot",
          xlab="thinned iteration number",ylab=expression(bold(U)))
  matplot(t(rot.loadings.store[sample(1:p,ceiling(.3*p)),sample(1:q,1),burn:nthin]), type="l",main="Loadings trace plot",
          xlab="thinned iteration number",ylab=expression(bold(W)))
  matplot(t(rot.beta.store[sample(1:q,1),,burn:nthin]), type="l",main="Regression coefficients trace plot",
          xlab="thinned iteration number",ylab=expression(bold(beta)))
  matplot(t(res2$sig.store[sample(1:p,ceiling(.3*p)),burn:nthin]), type="l",main="Covariance trace plot",
          xlab="thinned iteration number",ylab=expression(bold(sigma^2)))
  matplot(t(rot.s2.store[,burn:nthin]), type="l",main="Latent covariance trace plot",
          xlab="thinned iteration number",ylab=expression(Sigma))
  
  ##############################################
  #### Loadings and regression coefficient plots
  dload<-apply(rot.loadings.store[,,burn:nthin],c(1,2),median)/apply(rot.loadings.store[,,burn:nthin],c(1,2),sd)
  plot(dload[,1:2], xlab="factor:1", ylab="factor:2", type="n",main="Standardized loadings plot")
  text(dload[,1:2], colnames(Y), cex=1)
  abline(h=c(-2,2), col="red", lty=2)
  abline(v=c(-2,2), col="red", lty=2)
  dbeta<-apply(rot.beta.store[,,burn:nthin],c(1,2), median)/apply(rot.beta.store[,,burn:nthin],c(1,2), sd)
  plot(dbeta[1:2,],xlim=c(min(dbeta[1,])-.5,max(dbeta[1,])+.5),ylim=c(min(dbeta[2,])-.5,max(dbeta[2,])+.5), 
       xlab="factor:1", ylab="factor:2", type="n",main="Standardized effects of covariates")
  text(dbeta[1,],dbeta[2,], colnames(X), cex=1.2)
  abline(h=c(-2,2), col="red", lty=2)
  abline(v=c(-2,2), col="red", lty=2)

  #######################################
  #### (111) K means and Model based clustering
  scores<-t(apply(rot.scores.store[,,burn:nthin,drop=FALSE],c(1,2), mean))
  cl<-kmeans(scores, q+1, nstart = 100)
  plot(scores[,1:2], col = cl$cluster,xlab="factor:1", ylab="factor:2",main="K means clustering")
  yBIC <- mclustBIC(scores,G=gmin:gmax)
  #plot(yBIC,main="Model based clustering",xlab="number of subtypes")
  rs1<-summary(yBIC,scores)
  tau<-rs1$classification
  plot(scores[,1:2],col=tau,xlab="factor:1", ylab="factor:2",main="Model based clustering")
  mb_km_t<-table(tau,cl$cluster)
  colnames(mb_km_t)<-paste("kmeans",1:length(table(cl$cluster)),sep="")
  rownames(mb_km_t)<-paste("mbc",1:length(table(tau)),sep="")
  mb_km_t
  
  ############
  #### Output
  list(b.store=rot.beta.store,u.store=rot.scores.store,w.store=rot.loadings.store,s.store=res2$sig.store, 
       s2.store=rot.s2.store,KMeansClustering=cl$cluster,MBClustering=tau,mcmc=mcmc,burn=burn,thinstep=thinstep,qBIC=res1$BIC)
}
