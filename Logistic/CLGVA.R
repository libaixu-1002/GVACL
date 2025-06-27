#library(statmod)
#source("CalculateB.Rs")
###############################################################################

ReGVA.CLCR.f <- function(vtheta,family,vy,mX,gh){
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  ###########################################################################
  vbetar <- vtheta[1]
  vbetac <- vtheta[2]
  vbeta1 <- vtheta[3:(p+1)]
  tu <- vtheta[p+2]
  tv <- vtheta[p+3]
  vym <- vy
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta1[i]
  }
  
  vyxr <- matrix(vbetar,m,n) + vyx0
  vyxc <- matrix(vbetac,m,n) + vyx0
  #dim(vyx0)
  ###########################################################################
  elbo <- 0
  for (i in 1:m) {
    muui <- as.numeric(vtheta[p+2*i+2])
    tui  <- as.numeric(vtheta[p+2*i+3])
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    vB0 <- B0rc.fun(family,vmu,vsigma2,gh)
    
    elbo <- elbo + sum(vym[i,]*vmu) - sum(vB0) - 0.5*tu - 0.5*(muui^2 + exp(tui))*exp(-tu)  + 0.5*tui + 0.5
  }
  
  for(j in 1:n){
    muvj <- as.numeric(vtheta[p+2*m+2*j+2])
    tvj  <- as.numeric(vtheta[p+2*m+2*j+3])
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    vB0 <- B0rc.fun(family,vmu,vsigma2,gh)
    
    elbo <- elbo + sum(vym[,j]*vmu) - sum(vB0) - 0.5*tv - 0.5*(muvj^2 + exp(tvj))*exp(-tv) + 0.5*tvj + 0.5
  }
  
  return(elbo)
}

###############################################################################

ReGVA.CLCR.vg <- function(vtheta,family,vy,mX,gh) 
{ 
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbetar <- vtheta[1]
  vbetac <- vtheta[2]
  vbeta1 <- vtheta[3:(p+1)]
  tu <- vtheta[p+2]
  tv <- vtheta[p+3]
  
  vym <- vy
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta1[i]
  }
  #for(i in 1:m){
  #  for(j in 1:n){
  #    vyx0[i,j] <- sum(mX[i,j,]*vbeta1)
  #  }
  #}
  vyxr <- matrix(vbetar,m,n) + vyx0
  vyxc <- matrix(vbetac,m,n) + vyx0
  ###########################################################################
  
  ELBO_betar <- 0
  ELBO_betac <- 0
  ELBO_beta1 <- rep(0,p-1)
  ELBO_tu <- 0
  ELBO_tv <- 0
  
  lvg.uxi <- list()
  
  for(i in 1:m){
    muui <- as.numeric(vtheta[p+2*i+2])
    tui  <- as.numeric(vtheta[p+2*i+3])
    
    ELBO_muui <- 0
    ELBO_tui  <- 0
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    #ELBO_betar <- ELBO_betar + sum(vym[i,]-res$vB1)
    #ELBO_beta1 <- ELBO_beta1 + sum((vym[i,]-res$vB1)*mX[i,])
    #ELBO_muui <-  sum(vym[i,]-res$vB1) - muui*exp(-tu)
    #ELBO_tui  <- -0.5*exp(tui)*sum(res$vB2) - 0.5*exp(tui)*exp(-tu) + 0.5
    
    
    midsum1 <- sum(vym[i,]-res$vB1)
    ELBO_betar <- ELBO_betar + midsum1
    
    ##midtemp <- as.matrix((vym[i,]-res$vB1)*mX[i,,],n,p-1)
    ELBO_beta1 <- ELBO_beta1 +   c(t(mX[i,,])%*%(vym[i,]-res$vB1))
    ##apply(midtemp,2,sum)
    ELBO_muui <-  midsum1 - muui*exp(-tu)
    ELBO_tui  <- -0.5*exp(tui)*sum(res$vB2) - 0.5*exp(tui)*exp(-tu) + 0.5
    
    ELBO_tu <- ELBO_tu - 0.5 + 0.5*(muui^2 + exp(tui))*exp(-tu)
    lvg.uxi[[i]] <- c(ELBO_muui,ELBO_tui)
  }
  
  
  lvg.vxj <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+2]
    tvj  <- vtheta[p+2*m+2*j+3]
    
    ELBO_muvj <- 0
    ELBO_tvj  <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    midsum2 <- sum(vym[,j] - res$vB1)
    ELBO_betac <- ELBO_betac + midsum2
    ELBO_beta1 <- ELBO_beta1 + c(t(mX[,j,])%*%(vym[,j]-res$vB1))
    ELBO_muvj <-  midsum2 - muvj*exp(-tv)
    ELBO_tvj <-  -0.5*exp(tvj)*sum(res$vB2) - 0.5*exp(tvj)*exp(-tv) + 0.5
    
    ELBO_tv <- ELBO_tv - 0.5 + 0.5*(muvj^2 + exp(tvj))*exp(-tv)
    lvg.vxj[[j]] <- c(ELBO_muvj,ELBO_tvj)
  }
  
  vg <- c()
  vg.theta <- c(ELBO_betar,ELBO_betac,ELBO_beta1,ELBO_tu,ELBO_tv)
  vg <- c(vg.theta)
  for(i in 1:m){
    vg <-c(vg,lvg.uxi[[i]])
  }
  for(j in 1:n){
    vg <-c(vg,lvg.vxj[[j]])
  }
  
  res.vg <- list(vg.0=c(vg.theta),lvg.i=lvg.uxi,lvg.j= lvg.vxj,vg=vg)
  
  return(res.vg)
}

###############################################################################
ReGVA.CLCR.vgANDmH <- function(vtheta,family,vy,mX,gh){ 
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  ###########################################################################
  vbetar <- vtheta[1]
  vbetac <- vtheta[2]
  vbeta1 <- vtheta[3:(p+1)]
  tu <- vtheta[p+2]
  tv <- vtheta[p+3]
  
  vym <- vy
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta1[i]
  }
  
  vyxr <- matrix(vbetar,m,n) + vyx0
  vyxc <- matrix(vbetac,m,n) + vyx0
  ###########################################################################
  
  ELBO_betar <- 0
  ELBO_betac <- 0
  ELBO_beta1 <- rep(0,p-1)
  ELBO_tu <- 0
  ELBO_tv <- 0
  
  ###########################################################################
  
  res.vg <- ReGVA.CLCR.vg(vtheta,family,vy,mX,gh) 
  
  ###########################################################################
  mH.theta <- matrix(0,p+3,p+3)
  ELBO_betabeta <- matrix(0,p+1,p+1)
  ELBO_tutu <- 0
  ELBO_tvtv <- 0
  
  for(i in 1:m){
    muui <- vtheta[p+2*i+2]
    tui  <- vtheta[p+2*i+3]
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ELBO_betabeta[1,1] <- ELBO_betabeta[1,1] - sum(res$vB2)
    #midtemp <- as.matrix(mX[i,,]*res$vB2,n,p-1)
    ELBO_betabeta[1,3:(p+1)] <- ELBO_betabeta[1,3:(p+1)] - c(t(mX[i,,])%*%res$vB2)
    #apply(midtemp,2,sum) 
    ELBO_betabeta[3:(p+1),3:(p+1)] <- ELBO_betabeta[3:(p+1),3:(p+1)] - t(mX[i,,])%*%(mX[i,,]*res$vB2) 
    
    ELBO_tutu <- ELBO_tutu - 0.5*(muui^2 + exp(tui))*exp(-tu)
  }
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+2]
    tvj  <- vtheta[p+2*m+2*j+3]
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ELBO_betabeta[2,2] <- ELBO_betabeta[2,2] - sum(res$vB2)
    ##midtemp <- as.matrix(mX[,j,]*res$vB2,m,p-1)
    ELBO_betabeta[2,3:(p+1)] <- ELBO_betabeta[2,3:(p+1)] - c(t(mX[,j,])%*%res$vB2)
    ##apply(midtemp,2,sum)
    ELBO_betabeta[3:(p+1),3:(p+1)] <- ELBO_betabeta[3:(p+1),3:(p+1)] - t(mX[,j,])%*%(mX[,j,]*res$vB2)
    
    ELBO_tvtv <- ELBO_tvtv - 0.5*(muvj^2 + exp(tvj))*exp(-tv)
  }
  
  ELBO_betabeta[3:(p+1),1] <- ELBO_betabeta[1,3:(p+1)]
  ELBO_betabeta[3:(p+1),2] <- ELBO_betabeta[2,3:(p+1)]
  
  mH.theta[1:(p+1),1:(p+1)] <- ELBO_betabeta
  mH.theta[p+2,p+2] <- ELBO_tutu
  mH.theta[p+3,p+3] <- ELBO_tvtv
  
  ##################################################
  mH.vivi <- matrix(0,2,2)
  mH.wjwj <- matrix(0,2,2)
  
  lvxi.uu <- list()
  
  for(i in 1:m){
    
    muui <- vtheta[p+2*i+2]
    tui  <- vtheta[p+2*i+3]
    
    ELBO_muuimuui <- 0
    ELBO_muuitui  <- 0
    ELBO_tuitui   <- 0
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    #midsum3 <- sum(res$vB2)
    ELBO_muuimuui <- -sum(res$vB2) - exp(-tu)
    ELBO_muuitui  <- -sum(res$vB3)*0.5*exp(tui)
    ELBO_tuitui   <- -sum(res$vB4)*0.25*exp(2*tui) - sum(res$vB2)*0.5*exp(tui) - exp(-tu)*exp(tui)*0.5 
    
    mH.uiui <- matrix(c(ELBO_muuimuui,ELBO_muuitui,ELBO_muuitui,ELBO_tuitui),2,2)
    lvxi.uu[[i]] <- mH.uiui
  }
  
  lwxi.vv <- list()
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+2]
    tvj  <- vtheta[p+2*m+2*j+3]
    
    ELBO_muvjmuvj <- 0
    ELBO_muvjtvj  <- 0
    ELBO_tvjtvj   <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ##midsum4 <- sum(res$vB2)
    ELBO_muvjmuvj <- -sum(res$vB2) - exp(-tv) 
    ELBO_muvjtvj  <- -sum(res$vB3)*0.5*exp(tvj)
    ELBO_tvjtvj   <- -sum(res$vB4)*0.25*exp(2*tvj) - sum(res$vB2)*0.5*exp(tvj) - exp(-tv)*exp(tvj)*0.5
    
    mH.vjvj <- matrix(c(ELBO_muvjmuvj,ELBO_muvjtvj,ELBO_muvjtvj,ELBO_tvjtvj),2,2)
    
    lwxi.vv[[j]] <- mH.vjvj
  }
  
  lwxi.uivj <- list()
  lxi.uv <- list()
  mH.xiuv <- matrix(0,2*(m+n),2*(m+n))
  mH.uu <- matrix(0,2*m,2*m)
  mH.vv <- matrix(0,2*n,2*n)
  for(i in 1:m){
    mH.uu[(2*i-1):(2*i),(2*i-1):(2*i)] <- lvxi.uu[[i]]
  }
  for(j in 1:n){
    mH.vv[(2*j-1):(2*j),(2*j-1):(2*j)] <- lwxi.vv[[j]]
  }
  
  mH.xiuv[1:(2*m),1:(2*m)] <- mH.uu
  mH.xiuv[1:(2*n)+2*m,1:(2*n)+2*m] <- mH.vv
  
  luthetaxi <- list()
  for(i in 1:m){
    muui <- vtheta[p+2*i+2]
    tui  <- vtheta[p+2*i+3]
    
    ELBO_betarmuui <- 0
    ELBO_beta1muui <- rep(0,p-1)
    ELBO_tumuui <- 0
    
    ELBO_betartui <- 0
    ELBO_beta1tui <- rep(0,p-1)
    ELBO_tutui <- 0
    
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ELBO_betarmuui <- -sum(res$vB2)
    
    #midtemp <- as.matrix(mX[i,,]*res$vB2,n,p-1)
    ELBO_beta1muui <- -c(t(mX[i,,])%*%res$vB2)
    #apply(midtemp, 2, sum)
    ELBO_betartui  <- -sum(res$vB3)*0.5*exp(tui)
    
    #midtemp <- as.matrix(mX[i,,]*res$vB3,n,p-1)
    ELBO_beta1tui  <- -c(t(mX[i,,])%*%res$vB3*0.5*exp(tui))
    #apply(midtemp,2,sum)*0.5*exp(tui)
    
    ELBO_tumuui <- muui*exp(-tu)
    ELBO_tutui  <- 0.5*exp(tui)*exp(-tu)
    
    mH.vitheta  <- matrix(c(ELBO_betarmuui,0,ELBO_beta1muui,ELBO_tumuui,0,
                            ELBO_betartui ,0,ELBO_beta1tui ,ELBO_tutui ,0),p+3,2)
    
    luthetaxi[[i]] <- mH.vitheta
  }
  
  lvthetaxi <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+2]
    tvj  <- vtheta[p+2*m+2*j+3]
    
    ELBO_betacmuvj <- 0
    ELBO_beta1muvj <- rep(0,p-1)
    ELBO_tvmuvj <- 0
    
    ELBO_betactvj <- 0
    ELBO_beta1tvj <- rep(0,p-1)
    ELBO_tvtvj <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    
    ELBO_betacmuvj <- -sum(res$vB2)
    
    #midtemp <- as.matrix(mX[,j,]*res$vB2,m,p-1)
    ELBO_beta1muvj <- -c(t(mX[,j,])%*%res$vB2)
    #apply(midtemp, 2, sum)
    ELBO_betactvj  <- -sum(res$vB3)*0.5*exp(tvj)
    ELBO_beta1tvj  <- -c(t(mX[,j,])%*%res$vB3)*0.5*exp(tvj)
    
    ELBO_tvmuvj <- muvj*exp(-tv)
    ELBO_tvtvj  <- 0.5*exp(tvj)*exp(-tv)
    
    mH.wjtheta  <- matrix(c(0,ELBO_betacmuvj,ELBO_beta1muvj,0,ELBO_tvmuvj,
                            0,ELBO_betactvj ,ELBO_beta1tvj ,0,ELBO_tvtvj),p+3,2)
    
    lvthetaxi[[j]] <- mH.wjtheta
  }
  
  Hthetaxi <- matrix(0,p+3,2*m + 2*n)
  Hthetau <- matrix(0,p+3,2*m)
  Hthetav <- matrix(0,p+3,2*n)
  
  for(i in 1:m){
    Hthetau[,(2*i-1):(2*i)] <- luthetaxi[[i]] 
  }
  for(j in 1:n){
    Hthetav[,(2*j-1):(2*j)] <- lvthetaxi[[j]] 
  }
  Hthetaxi[,1:(2*m)] <- Hthetau
  Hthetaxi[,1:(2*n)+2*m] <- Hthetav
  
  res.mH <- list(mH.00=mH.theta,lmH.0i=Hthetaxi,lmH.ii=mH.xiuv)
  res.vg <- ReGVA.CLCR.vg(vtheta,family,vy,mX,gh) 
  
  return(list(res.vg=res.vg,res.mH=res.mH))
}

ReGVA.CLCR.NRDIR <- function(vtheta,family,vy,mX,gh) 
{
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbetar <- vtheta[1]
  vbetac <- vtheta[2]
  vbeta1 <- vtheta[3:(p+1)]
  tu <- vtheta[p+2]
  tv <- vtheta[p+3]
  
  vym <- vy
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta1[i]
  }
  
  vyxr <- matrix(vbetar,m,n) + vyx0
  vyxc <- matrix(vbetac,m,n) + vyx0
  ###########################################################################
  
  # Calculate 1d integrals and some auxiliary information
  res <- ReGVA.CLCR.vgANDmH(vtheta,family,vy,mX,gh) 
  res.vg <- res$res.vg
  res.mH <- res$res.mH
  vg.til <- res.vg$vg.0
  suv <- res.vg$vg[-c(1:(p+3))]
  
  mH.til <- res.mH$mH.00
  
  mH.thetavw <- res.mH$lmH.0i
  mH.xivw <- res.mH$lmH.ii
  
  mHxivw.inv <- solve(mH.xivw)
  mH.til <- mH.til - mH.thetavw%*%mHxivw.inv%*%t(mH.thetavw)
  vg.til <- vg.til - mH.thetavw%*%mHxivw.inv%*%suv
  
  sI <- solve(mH.til)
  vs <- sI%*%vg.til
  
  dvtheta <- vs
  dvtheta <- c(dvtheta,solve(res.mH$lmH.ii,suv - t(mH.thetavw)%*%vs))
  
  return(list(dvtheta=-dvtheta,vg=res.vg$vg,sI=-sI))
}


###############################################################################
ReGVA.CLCR.FIT <- function(vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id) 
{
  MAXITER  <- 100
  MAXJTER  <- 20
  EPS.TERM.PAR  <- 5.0E-4
  EPS.TERM.GRAD <- 1.0E-9
  
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  
  tuv <- log(vSigma)
  uxi <- c()
  for(i in 1:m){
    tui <- log(luLambda[[i]])
    uxi <- c(uxi,lumu[[i]],tui)
  }
  vxi <- c()
  for(j in 1:n){
    tvj <- log(lvLambda[[j]])
    vxi <- c(vxi,lvmu[[j]],tvj)
  }
  vbeta <- rep(0,p+1)
  vtheta <- c(vbeta,tuv,uxi,vxi)
  sl <- ReGVA.CLCR.f(vtheta,family,vy,mX,gh) 
  NI <- 0
  cat("Initial variational log-compositelikelihood = ",sl)
  cat("\n\n")
  
  cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
  
  for (ITER in 1:MAXITER) {
    vtheta.old <- vtheta
    
    # Calculate Newton-Raphson Direction 
    res <- ReGVA.CLCR.NRDIR(vtheta,family,vy,mX,gh)    
    
    stepsize <- 1                
    vtheta.temp <- vtheta + stepsize*res$dvtheta
    sl.temp <- ReGVA.CLCR.f(vtheta.temp,family,vy,mX,gh) 
    if (is.nan(sl.temp)) { 
      DOGRIDSEARCH <- TRUE
    } else {            
      if (sl.temp>sl) { DOGRIDSEARCH <- FALSE } 
      else { DOGRIDSEARCH <- TRUE }
    }         
    
    if (DOGRIDSEARCH) {
      NG <- 15
      vsl <- c()
      
      # Allow for small negative step sizes 
      # (just in case Hessian is negitive definite)
      vstepsize <- seq(-0.25,1,,NG)[-NG]  
      
      for (JTER in 1:length(vstepsize)) {
        vtheta.temp <- vtheta + vstepsize[JTER]*res$dvtheta
        vsl[JTER] <- ReGVA.CLCR.f(vtheta.temp,family,vy,mX,gh) 
      }                   
      sl <- max(vsl)
      stepsize <- vstepsize[which.max(vsl)]
      vtheta <- vtheta + stepsize*res$dvtheta
    } else {
      sl <- sl.temp
      vtheta <- vtheta.temp           
    }
    
    
    # Check termination condtions
    err.par <- max(abs(vtheta-vtheta.old))/max(abs(vtheta))
    err.vg  <- sum(res$vg^2)/length(vtheta)
    
    cat(sprintf("%4d    %E    %+E    %E  %E\n",ITER,sl,stepsize,err.vg,err.par))
    NI <- ITER
    if (err.vg < EPS.TERM.GRAD) { break; }
    if (err.par < EPS.TERM.PAR) { break; }
    
    if (exp(vtheta[p+2])<EPS.TERM.PAR) { break; } # Break if sigma2 is getting too small
    if (exp(vtheta[p+3])<EPS.TERM.PAR) { break; }
  }
  
  vbeta <- vtheta[1:(p+1)]
  sigmauv <- vtheta[p+2:3]
  
  vmuu   <- vtheta[p+3+2*(1:m)-1]
  vzetau <- vtheta[p+3+2*(1:m)]
  
  vmuv<- vtheta[p+3+ 2*m + 2*(1:n)-1]
  vzetav <- vtheta[p+3+ 2*m + 2*(1:n)]
  
  return(list(NI = NI,vbeta=vbeta,sigma2=exp(sigmauv),muu=vmuu,lambdau=exp(vzetau),muv=vmuv,lambdav=exp(vzetav),vbeta.serr=sqrt(abs(diag(res$sI))),sl=sl))
}

###############################################################################

