
#library(statmod)
#source("CalculateB.Rs")

###############################################################################

ReGVA.CLCR.f <- function(alpha,vtheta,family,vy,mX,gh){
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
  elbo <- 0
  #elbo_1 <- 0
  #elbo_2 <- 0
  #elbo_3 <- 0
  #elbo_4 <- 0
  #elbo_5 <- 0
  #elbo_6 <- 0
  
  #vmu.X  <- mX%*%vbeta
  for (i in 1:m) {
    muui <- as.numeric(vtheta[p+(i-1)*2+4])
    tui  <- as.numeric(vtheta[p+(i-1)*2+5])
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    vB0 <- B0rc.fun(family,vmu,vsigma2,gh)
    
    cij <- n*alpha*log(alpha) + (alpha-1)*sum(log(vym[i,])) - n*lgamma(alpha)
    
    elbo <- elbo - alpha*sum(vym[i,]*vB0)- alpha*sum(vmu) + cij
    elbo <- elbo - 0.5*log(2*pi) - 0.5*tu - 0.5*(muui^2 + exp(tui))*exp(-tu)
    elbo <- elbo + 0.5*log(2*pi) + 0.5*tui + 0.5
  }
  
  for(j in 1:n){
    muvj <- as.numeric(vtheta[p+2*m+(j-1)*2+4])
    tvj  <- as.numeric(vtheta[p+2*m+(j-1)*2+5])
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    vB0 <- B0rc.fun(family,vmu,vsigma2,gh)
    
    cij <- m*alpha*log(alpha) + (alpha-1)*sum(log(vym[,j])) - m*lgamma(alpha)
    
    elbo <- elbo - alpha*sum(vym[,j]*vB0) - alpha*sum(vmu) + cij
    elbo <- elbo - 0.5*log(2*pi) - 0.5*tv - 0.5*(muvj^2 + exp(tvj))*exp(-tv)
    elbo <- elbo + 0.5*log(2*pi) + 0.5*tvj + 0.5
  }
  
  #elbo <- elbo_1 + elbo_2 + elbo_3 + elbo_4 + elbo_5 + elbo_6
  
  return(elbo)
}

###############################################################################

ReGVA.CLCR.vg <- function(alpha,vtheta,family,vy,mX,gh) 
{ 
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
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
  
  lvg.uxi <- list()
  
  for(i in 1:m){
    muui <- as.numeric(vtheta[p+(i-1)*2+4])
    tui  <- as.numeric(vtheta[p+(i-1)*2+5])
    
    ELBO_muui <- 0
    ELBO_tui  <- 0
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    #midsum1 <- alpha*sum(vym[i,]*res$vB1-1)
    
    ELBO_betar <- ELBO_betar + alpha*sum(vym[i,]*res$vB1-1)
    
    #midtemp <- as.matrix((vym[i,]*res$vB1-1)*mX[i,,],n,p-1)
    #apply(midtemp, 2, sum)
    ELBO_beta1 <- ELBO_beta1 + alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB1-1))
    #sum((vym[i,]*res$vB1-1)*mX[i,,])
    
    ELBO_muui  <-  alpha*sum(vym[i,]*res$vB1-1) - muui*exp(-tu)
    ELBO_tui   <- -0.5*exp(tui)*alpha*sum(vym[i,]*res$vB2) - 0.5*exp(tui)*exp(-tu) + 0.5
    
    ELBO_tu <- ELBO_tu - 0.5 + 0.5*(muui^2 + exp(tui))*exp(-tu)
    lvg.uxi[[i]] <- c(ELBO_muui,ELBO_tui)
  }
  
  
  lvg.vxj <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+(j-1)*2+4]
    tvj  <- vtheta[p+2*m+(j-1)*2+5]
    
    ELBO_muvj <- 0
    ELBO_tvj  <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    #midsum2 <- alpha*sum(vym[,j]*res$vB1-1) 
    ELBO_betac <- ELBO_betac + alpha*sum(vym[,j]*res$vB1-1) 
    
    #midtemp <-  as.matrix((vym[,j]*res$vB1-1)*mX[,j,],m,p-1)
    ELBO_beta1 <- ELBO_beta1 + alpha*c(t(mX[,j,])%*%(vym[,j]*res$vB1-1))
    #apply(midtemp,2,sum)
    
    ELBO_muvj <-  alpha*sum(vym[,j]*res$vB1-1) - muvj*exp(-tv)
    ELBO_tvj  <-  -0.5*exp(tvj)*alpha*sum(vym[,j]*res$vB2) - 0.5*exp(tvj)*exp(-tv) + 0.5
    
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
ReGVA.CLCR.vgANDmH <- function(alpha,vtheta,family,vy,mX,gh){ 
  
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
  # EPS <- 1.0E-15
  ###########################################################################
  
  res.vg <- ReGVA.CLCR.vg(alpha,vtheta,family,vy,mX,gh) 
  
  ###########################################################################
  mH.theta <- matrix(0,p+3,p+3)
  ELBO_betabeta <- matrix(0,p+1,p+1)
  ELBO_tutu <- 0
  ELBO_tvtv <- 0
  
  for(i in 1:m){
    muui <- vtheta[p+(i-1)*2+4]
    tui  <- vtheta[p+(i-1)*2+5]
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ELBO_betabeta[1,1] <- ELBO_betabeta[1,1] - alpha*sum(vym[i,]*res$vB2) 
    #midtemp <- as.matrix(mX[i,,]*vym[i,]*res$vB2,n,p-1)
    ELBO_betabeta[1,3:(p+1)] <- ELBO_betabeta[1,3:(p+1)] - alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB2))
    #apply(midtemp, 2, sum)
    #alpha*sum(vym[i,]*mX[i,]*res$vB2) 
    ELBO_betabeta[3:(p+1),3:(p+1)] <- ELBO_betabeta[3:(p+1),3:(p+1)] - alpha*(t(mX[i,,])%*%(mX[i,,]*vym[i,]*res$vB2))
    #alpha*sum(mX[i,]^{2}*vym[i,]*res$vB2) 
    
    ELBO_tutu <- ELBO_tutu - 0.5*(muui^2 + exp(tui))*exp(-tu)
  }
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+(j-1)*2+4]
    tvj  <- vtheta[p+2*m+(j-1)*2+5]
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    ELBO_betabeta[2,2] <- ELBO_betabeta[2,2] - alpha*sum(vym[,j]*res$vB2)
    #midtemp <- as.matrix(mX[,j,]*vym[,j]*res$vB2)
    ELBO_betabeta[2,3:(p+1)] <- ELBO_betabeta[2,3:(p+1)] - alpha*c(t(mX[,j,])%*%(vym[,j]*res$vB2))
    #apply(midtemp, 2, sum)
    #alpha*sum(mX[,j]*vym[,j]*res$vB2)
    ELBO_betabeta[3:(p+1),3:(p+1)] <- ELBO_betabeta[3:(p+1),3:(p+1)] - alpha*(t(mX[,j,])%*%(mX[,j,]*vym[,j]*res$vB2))
      
    #alpha*sum(mX[,j]^{2}*vym[,j]*res$vB2)
    
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
    
    muui <- vtheta[p+(i-1)*2+4]
    tui  <- vtheta[p+(i-1)*2+5]
    
    ELBO_muuimuui <- 0
    ELBO_muuitui  <- 0
    ELBO_tuitui   <- 0
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    midsum3 <- alpha*sum(vym[i,]*res$vB2)
    ELBO_muuimuui <- -midsum3 - exp(-tu)
    ELBO_muuitui  <-  midsum3*0.5*exp(tui)
    ELBO_tuitui   <- -midsum3*0.25*exp(2*tui) - midsum3*0.5*exp(tui) - exp(-tu)*exp(tui)*0.5 
    
    mH.uiui <- matrix(c(ELBO_muuimuui,ELBO_muuitui,ELBO_muuitui,ELBO_tuitui),2,2)
    lvxi.uu[[i]] <- mH.uiui
  }
  
  lwxi.vv <- list()
  for(j in 1:n){
    muvj <- vtheta[p+2*m+(j-1)*2+4]
    tvj  <- vtheta[p+2*m+(j-1)*2+5]
    
    ELBO_muvjmuvj <- 0
    ELBO_muvjtvj  <- 0
    ELBO_tvjtvj   <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    midsum4 <- alpha*sum(vym[,j]*res$vB2)
    ELBO_muvjmuvj <- -midsum4 - exp(-tv) 
    ELBO_muvjtvj  <-  midsum4*0.5*exp(tvj)
    ELBO_tvjtvj   <- -midsum4*0.25*exp(2*tvj) - midsum4*0.5*exp(tvj) - exp(-tv)*exp(tvj)*0.5
    
    mH.vjvj <- matrix(c(ELBO_muvjmuvj,ELBO_muvjtvj,ELBO_muvjtvj,ELBO_tvjtvj),2,2)
    
    lwxi.vv[[j]] <- mH.vjvj
  }
  
  lwxi.uivj <- list()
  lxi.uv <- list()
  mH.xiuv <- matrix(0,2*(m+n),2*(m+n))
  mH.uu <- matrix(0,2*m,2*m)
  mH.vv <- matrix(0,2*n,2*n)
  #mH.uv <- matrix(0,2*m,2*ni)
  for(i in 1:m){
    mH.uu[(2*i-1):(2*i),(2*i-1):(2*i)] <- lvxi.uu[[i]]
  }
  for(j in 1:n){
    mH.vv[(2*j-1):(2*j),(2*j-1):(2*j)] <- lwxi.vv[[j]]
  }
  
  mH.xiuv[1:(2*m),1:(2*m)] <- mH.uu
  mH.xiuv[1:(2*n)+2*m,1:(2*n)+2*m] <- mH.vv
  #mH.xivw[1:(2*m),1:(2*ni)+2*m] <- mH.vw
  #mH.xivw[1:(2*ni)+2*m,1:(2*m)] <- t(mH.vw)
  
  luthetaxi <- list()
  for(i in 1:m){
    muui <- vtheta[p+(i-1)*2+4]
    tui  <- vtheta[p+(i-1)*2+5]
    
    ELBO_betarmuui <- 0
    ELBO_beta1muui <- rep(0,p-1)
    ELBO_tumuui <- 0
    
    ELBO_betartui <- 0
    ELBO_beta1tui <- rep(0,p-1)
    ELBO_tutui <- 0
    
    
    vmu <- vyxr[i,] + muui
    vsigma2 <- rep(exp(tui),n)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    
    ELBO_betarmuui <- -alpha*sum(vym[i,]*res$vB2)
    
    #midtemp <- as.matrix(vym[i,]*mX[i,,]*res$vB2,n,p-1)
    ELBO_beta1muui <- -alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB2))
    #apply(midtemp, 2, sum)
    #sum(vym[i,]*mX[i,]*res$vB2)
   
    ELBO_betartui  <-  alpha*sum(vym[i,]*res$vB3)*0.5*exp(tui)
    
    #midtemp <- as.matrix(vym[i,]*mX[i,,]*res$vB3,n,p-1)
    ELBO_beta1tui  <-  alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB3))*0.5*exp(tui)
    #apply(midtemp, 2, sum)*0.5*exp(tui)
    
    
    ELBO_tumuui <- muui*exp(-tu)
    ELBO_tutui  <- 0.5*exp(tui)*exp(-tu)
    
    mH.vitheta  <- matrix(c(ELBO_betarmuui,0,ELBO_beta1muui,ELBO_tumuui,0,ELBO_betartui ,0,ELBO_beta1tui ,ELBO_tutui ,0),p+3,2)
    
    luthetaxi[[i]] <- mH.vitheta
  }
  
  lvthetaxi <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+(j-1)*2+4]
    tvj  <- vtheta[p+2*m+(j-1)*2+5]
    
    ELBO_betacmuvj <- 0
    ELBO_beta1muvj <- rep(0,p-1)
    ELBO_tvmuvj <- 0
    
    ELBO_betactvj <- 0
    ELBO_beta1tvj <- rep(0,p-1)
    ELBO_tvtvj <- 0
    
    vmu <- vyxc[,j] + muvj
    vsigma2 <- rep(exp(tvj),m)
    
    res <- B1234rc.fun(family,vmu,vsigma2,gh)
    
    
    ELBO_betacmuvj <- -alpha*sum(vym[,j]*res$vB2)
    #midtemp <- as.matrix(vym[,j]*mX[,j,]*res$vB2,m,p-1)
    ELBO_beta1muvj <- -alpha*c(t(mX[,j,])%*%(vym[,j]*res$vB2))
    #apply(midtemp,2,sum)
    #sum(vym[,j]*mX[,j]*res$vB2)
    ELBO_betactvj  <-  alpha*sum(vym[,j]*res$vB3)*0.5*exp(tvj)
    ELBO_beta1tvj  <-  alpha*c(t(mX[,j,])%*%(vym[,j]*res$vB3))*0.5*exp(tvj)
    #alpha*sum(vym[,j]*mX[,j]*res$vB3)*0.5*exp(tvj)
    #apply(midtemp, 2, sum)
    
    
    ELBO_tvmuvj <- muvj*exp(-tv)
    ELBO_tvtvj <- 0.5*exp(tvj)*exp(-tv)
    
    mH.wjtheta  <- matrix(c(0,ELBO_betacmuvj,ELBO_beta1muvj,0,ELBO_tvmuvj,0,ELBO_betactvj ,ELBO_beta1tvj ,0,ELBO_tvtvj),p+3,2)
    
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
  res.vg <- ReGVA.CLCR.vg(alpha,vtheta,family,vy,mX,gh) 
  
  return(list(res.vg=res.vg,res.mH=res.mH))
}


ReGVA.CLCR.NRDIR <- function(alpha,vtheta,family,vy,mX,gh) 
{
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
  
  #EPS <- 1.0E-15
  ###########################################################################
  
  # Calculate 1d integrals and some auxiliary information
  res <- ReGVA.CLCR.vgANDmH(alpha,vtheta,family,vy,mX,gh) 
  res.vg <- res$res.vg
  res.mH <- res$res.mH
  vg.til <- res.vg$vg.0
  suv <- res.vg$vg[-c(1:(p+3))]
  
  mH.til <- res.mH$mH.00
  
  mH.thetavw <- res.mH$lmH.0i
  mH.xivw <- res.mH$lmH.ii
  #dim(mH.xivw)
  #mH.xivw[1:10,1:10]
  #dim(res.mH$lmH.ii)
  
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
ReGVA.CLCR.FIT <- function(alpha,vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id) 
{
  MAXITER  <- 20
  MAXJTER  <- 20
  EPS.TERM.PAR  <- 1.0E-5
  EPS.TERM.GRAD <- 1.0E-9
  
  #family="POISSON"
  #family="LOGISTIC"
  
  if (family=="LOGISTIC") {
  # Generate Quadrature Points
    N <- 9
    gh.temp <- gauss.quad(N,kind="hermite")
    gh <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2))  
  }  
  
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  
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
  
  #vym <- vy
  #alpha <- 10
  sl <- ReGVA.CLCR.f(alpha,vtheta,family,vy,mX,gh) 
  
  cat("Initial variational log-compositelikelihood = ",sl)
  cat("\n\n")
  
  cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
  
  for (ITER in 1:MAXITER) {
    #ITER <- 1
    vtheta.old <- vtheta
    
    # Calculate Newton-Raphson Direction 
    res <- ReGVA.CLCR.NRDIR(alpha,vtheta,family,vy,mX,gh)    
    
    stepsize <- 1                
    vtheta.temp <- vtheta + stepsize*res$dvtheta
    sl.temp <- ReGVA.CLCR.f(alpha,vtheta.temp,family,vy,mX,gh) 
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
        vsl[JTER] <- ReGVA.CLCR.f(alpha,vtheta.temp,family,vy,mX,gh) 
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
    
    if (err.vg<EPS.TERM.GRAD) { break; }
    if (err.par<EPS.TERM.PAR) { break; }
    
    if (exp(vtheta[p+2])<5.0E-9) { break; } # Break if sigma2 is getting too small
    if (exp(vtheta[p+3])<5.0E-9) { break; }
  }
  
  vbeta <- vtheta[1:(p+1)]
  sigmauv <- vtheta[p+2:3]
  
  vmuu   <- vtheta[p+3+2*(1:m)-1]
  vzetau <- vtheta[p+3+2*(1:m)]
  
  vmuv<- vtheta[p+3 + 2*m + 2*(1:n)-1]
  vzetav <- vtheta[p+3 + 2*m + 2*(1:n)]
  
  return(list(vbeta=vbeta,sigma2=exp(sigmauv),muu=vmuu,lambdau=exp(vzetau),muv=vmuv,lambdav=exp(vzetav),vbeta.serr=sqrt(diag(res$sI)[1:(p+3)]),sl=sl))
}

###############################################################################

