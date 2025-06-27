
#library(statmod)
#source("CalculateB.Rs")

###############################################################################
###Gamma Distribution with fixed alpha###
GVA.UiVj.f <- function(alpha,vtheta,family,vy,mX,gh) 
{
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  
  vyx  <- matrix(vbeta[1],m,n) + vyx0
  
  elbo <- 0
  
  for (i in 1:m) {
    muui <- vtheta[p + 2*i + 1]
    tui  <- vtheta[p + 2*i + 2]
    
    muvj <- vtheta[p + 2*m + 1 + 2*(1:n)]
    tvj  <- vtheta[p + 2*m + 2 + 2*(1:n)]
    
    vmuuv <- vyx[i,] 
    vB0 <- B0.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    #length(vB0)
    ci <- n*alpha*log(alpha) + (alpha-1)*sum(log(vym[i,])) - n*log(lgamma(alpha))
    elbo <- elbo - alpha*sum(vym[i,]*vB0) - alpha*sum((vmuuv + muui + muvj)) + ci
    elbo <- elbo - 0.5*log(2*pi) - 0.5*tu - 0.5*(muui^2 + exp(tui))*exp(-tu)
    elbo <- elbo + 0.5*log(2*pi) + 0.5*tui + 0.5
  }
  
  
  muvj <- vtheta[p + 2*m + 2*(1:n) + 1]
  tvj  <- vtheta[p + 2*m + 2*(1:n) + 2]
  
  elbo <- elbo  - 0.5*n*tv - 0.5*sum(muvj^2 + exp(tvj))*exp(-tv)
  elbo <- elbo  + 0.5*sum(tvj) + 0.5*n
  
  #for(j in 1:ni){
  #  muvj <- vtheta[4 + 2*m + (j-1)*2 + 1]
  #  tvj  <- vtheta[4 + 2*m + (j-1)*2 + 2]
  #  
  #  elbo_4 <- elbo_4 - 0.5*log(2*pi) - 0.5*tv - 0.5*(muvj^2 + exp(tvj))*exp(-tv)
  #  elbo_5 <- elbo_5 + 0.5*log(2*pi) + 0.5*tvj + 0.5
  #}
  #elbo <- elbo_1 + elbo_2 + elbo_3 + elbo_4 + elbo_5
  
  return(elbo)
}

###############################################################################

GVA.UiVj.vg <- function(alpha,vtheta,family,vy,mX,gh) 
{ 
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  ###########################################################################
  vbeta <- vtheta[1:p]
  #vbeta0 <- vbeta[1]
  #vbeta1 <- vbeta[2]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]

  vym <- vy
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  
  vyx  <- matrix(vbeta[1],m,n) + vyx0
  
  ###########################################################################
  ELBO_beta <- rep(0,p)
  ELBO_tu <- 0
  ELBO_tv <- 0
  
  #for each i
  lvg.uxi <- list()
  for(i in 1:m){
    #i <- 1
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    ELBO_muui <- 0
    ELBO_tui  <- 0
    
    muvj <- vtheta[p + 2*m + 1 + 2*(1:n)]
    tvj  <- vtheta[p + 2*m + 2 + 2*(1:n)]
    
    vmuuv <- vyx[i,]
    
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_beta[1] <- ELBO_beta[1] + alpha*sum(vym[i,]*res$vB1-1)
    #midtemp <- (vym[i,]*res$vB1-1)*
    #  t(mX[i,,])%*%(vym[i,]*res$vB1-1)
    ELBO_beta[2:p] <- ELBO_beta[2:p] + alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB1-1))
    #  apply(midtemp,2,sum)
    ELBO_tu   <- ELBO_tu   - 0.5 + 0.5*(muui^2 + exp(tui))*exp(-tu)
    ELBO_muui <- ELBO_muui + alpha*sum(vym[i,]*res$vB1 - 1) - muui*exp(-tu)
    ELBO_tui  <- ELBO_tui - 0.5*alpha*exp(tui)*sum(vym[i,]*res$vB2) - 0.5*exp(tui)*exp(-tu) + 0.5
    
    lvg.uxi[[i]] <- c(ELBO_muui,ELBO_tui)
  }
  
  lvg.vxj <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+1]
    tvj  <- vtheta[p+2*m+2*j+2]
    
    ELBO_muvj <- 0
    ELBO_tvj  <- 0
    
    muui <- vtheta[p+2*(1:m)+1]
    tui  <- vtheta[p+2*(1:m)+2]
    
    vmuuv <- vyx[,j]
    
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    ELBO_muvj <- ELBO_muvj + alpha*sum(vym[,j]*res$vB1-1) - muvj*exp(-tv)
    ELBO_tvj <- ELBO_tvj - 0.5*alpha*exp(tvj)*sum(vym[,j]*res$vB2) - 0.5*exp(tvj)*exp(-tv) + 0.5
    ELBO_tv   <- ELBO_tv   - 0.5 + 0.5*(muvj^2 + exp(tvj))*exp(-tv)
    
    lvg.vxj[[j]] <- c(ELBO_muvj,ELBO_tvj)
  }
  
  vg <- c()
  vg.theta <- c(ELBO_beta,ELBO_tu,ELBO_tv)
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
GVA.UiVj.vgANDmH <- function(alpha,vtheta,family,vy,mX,gh){ 
  
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1
  
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  
  vyx  <- matrix(vbeta[1],m,n) + vyx0
  ###########################################################################
  res.vg <- GVA.UiVj.vg(alpha,vtheta,family,vy,mX,gh)
  ###########################################################################
  mH.theta <- matrix(0,p+2,p+2)
  ELBO_betabeta <- matrix(0,p,p)
  ELBO_tutu <- 0
  ELBO_tvtv <- 0
  
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    muvj <- vtheta[p + 2*m + 2*(1:n) + 1]
    tvj  <- vtheta[p + 2*m + 2*(1:n) + 2]
    
    vmuuv <- vyx[i,]
    
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_betabeta[1,1] <- ELBO_betabeta[1,1] - alpha*sum(vym[i,]*res$vB2)
    #midtemp <- mX[i,,]*(vym[i,]*res$vB2)
    #c(t(mX[i,,])%*%(vym[i,]*res$vB2))
    ELBO_betabeta[1,2:p] <- ELBO_betabeta[1,2:p] - alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB2))
    #apply(midtemp,2,sum)
    ###sum(mX[i,,]*(vym[i,]*res$vB2))
    #midtemp <- t(mX[i,,])%*%(mX[i,,]*(vym[i,]*res$vB2))
    ELBO_betabeta[2:p,2:p] <- ELBO_betabeta[2:p,2:p] - alpha*t(mX[i,,])%*%(mX[i,,]*(vym[i,]*res$vB2))
    ###alpha*sum(mX[i,]*mX[i,]*(vym[i,]*res$vB2))
    ELBO_betabeta[2:p,1] <- ELBO_betabeta[1,2:p]
    
    #ELBO_betabeta <- ELBO_betabeta - c(1,mX[i,,])
    
    ELBO_tutu <- ELBO_tutu - 0.5*(muui^2 + exp(tui))*exp(-tu)
  }
  
  
  muvj <- vtheta[p+2*m+2*(1:n)+1]
  tvj  <- vtheta[p+2*m+2*(1:n)+2]
  
  ELBO_tvtv <- ELBO_tvtv - 0.5*sum(muvj^2 + exp(tvj))*exp(-tv)
  
  mH.theta[1:p,1:p] <- ELBO_betabeta 
  #- diag(EPS,2)
  mH.theta[p+1,p+1] <- ELBO_tutu 
  mH.theta[p+2,p+2] <- ELBO_tvtv 
  
  ##################################################
  mH.uiui <- matrix(0,2,2)
  mH.vjvj <- matrix(0,2,2)
  
  lvxi.uu <- list()
  
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    ELBO_muuimuui <- 0
    ELBO_muuitui  <- 0
    ELBO_tuitui   <- 0
    
    muvj <- vtheta[p + 2*m + 2*(1:n) +1]
    tvj  <- vtheta[p + 2*m + 2*(1:n) +2]
    
    vmuuv <- vyx[i,]
    
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_muuimuui <- ELBO_muuimuui - alpha*sum(vym[i,]*res$vB2) - exp(-tu) 
    ELBO_muuitui  <- ELBO_muuitui  + alpha*sum(vym[i,]*res$vB3)*exp(tui)*0.5
    ELBO_tuitui   <- ELBO_tuitui   - alpha*sum(vym[i,]*res$vB4)*exp(2*tui)*0.25 - alpha*sum(vym[i,]*res$vB2)*exp(tui)*0.5 - exp(-tu)*exp(tui)*0.5 
    
    mH.uiui <- matrix(c(ELBO_muuimuui,ELBO_muuitui,ELBO_muuitui,ELBO_tuitui),2,2)
    lvxi.uu[[i]] <- mH.uiui
  }
  
  lwxi.vv <- list()
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+1]
    tvj  <- vtheta[p+2*m+2*j+2]
    
    ELBO_muvjmuvj <- 0
    ELBO_muvjtvj  <- 0
    ELBO_tvjtvj   <- 0
    
    muui <- vtheta[p + 2*(1:m) + 1]
    tui  <- vtheta[p + 2*(1:m) + 2]
    
    vmuuv <- vyx[,j]
    
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_muvjmuvj <- ELBO_muvjmuvj - alpha*sum(vym[,j]*res$vB2) - exp(-tv) 
    ELBO_muvjtvj  <- ELBO_muvjtvj  + alpha*sum(vym[,j]*res$vB3)*exp(tvj)*0.5
    ELBO_tvjtvj   <- ELBO_tvjtvj   - alpha*sum(vym[,j]*res$vB4)*exp(2*tvj)*0.25 - alpha*sum(vym[,j]*res$vB2)*exp(tvj)*0.5 - 0.5*exp(tvj)*exp(-tv)
    
    mH.vjvj <- matrix(c(ELBO_muvjmuvj,ELBO_muvjtvj,ELBO_muvjtvj,ELBO_tvjtvj),2,2) 
    
    lwxi.vv[[j]] <- mH.vjvj
    
  }
  
  lwxi.uivj <- list()
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    muvj <- vtheta[p + 2*m + 2*(1:n) + 1]
    tvj  <- vtheta[p + 2*m + 2*(1:n) + 2]
    
    vmuuv <- vyx[i,]
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    midsum1 <- alpha*vym[i,]*res$vB2
    ELBO_muuimuvj <- -midsum1
    ELBO_muuitvj  <- midsum1*0.5*exp(tvj)
    ELBO_tuimuvj  <- midsum1*0.5*exp(tui)
    ELBO_tuitvj   <- -midsum1*0.25*exp(tui+tvj)
    for(j in 1:n){
      index <- (i-1)*n + j
      lwxi.uivj[[index]] <- matrix(c(ELBO_muuimuvj[j],ELBO_tuimuvj[j],ELBO_muuitvj[j],ELBO_tuitvj[j]),2,2)
    }
  }
  
  mH.xiuv <- matrix(0,2*(m+n),2*(m+n))
  mH.uu <- matrix(0,2*m,2*m)
  mH.vv <- matrix(0,2*n,2*n)
  mH.uv <- matrix(0,2*m,2*n)
  
  for(i in 1:m){
    mH.uu[(2*i-1):(2*i),(2*i-1):(2*i)] <- lvxi.uu[[i]]
  }
  
  for(j in 1:n){
    mH.vv[(2*j-1):(2*j),(2*j-1):(2*j)] <- lwxi.vv[[j]]
  }
  
  for(i in 1:m){
    for(j in 1:n){
      index <- (i-1)*n + j
      mH.uv[(2*i-1):(2*i),(2*j-1):(2*j)] <- lwxi.uivj[[index]]
    }
  }
  
  mH.xiuv[1:(2*m),1:(2*m)] <- mH.uu
  mH.xiuv[1:(2*n)+2*m,1:(2*n)+2*m] <- mH.vv
  mH.xiuv[1:(2*m),(1:(2*n))+2*m] <- mH.uv
  mH.xiuv[(1:(2*n))+2*m,1:(2*m)] <- t(mH.uv)
  
  luthetaxi <- list()
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    ELBO_betamuui <- rep(0,p)
    ELBO_tumuui <- 0
    ELBO_tvmuui <- 0
    
    ELBO_betatui <- rep(0,p)
    ELBO_tutui <- 0
    ELBO_tvtui <- 0
    
    muvj <- vtheta[p+2*m+2*(1:n)+1]
    tvj  <- vtheta[p+2*m+2*(1:n)+2]
    
    vmuuv <- vyx[i,]
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_betamuui[1] <- ELBO_betamuui[1] - alpha*sum(vym[i,]*res$vB2)
    #midtemp <- mX[i,,]*(vym[i,]*res$vB2)
    #t(mX[i,,])%*%(vym[i,]*res$vB2)
    ELBO_betamuui[2:p] <- ELBO_betamuui[2] - alpha*c(t(mX[i,,])%*%(vym[i,]*res$vB2))
    ELBO_betatui[1]  <- ELBO_betatui[1]  + 0.5*alpha*exp(tui)*sum(vym[i,]*res$vB3)
    #midtemp <- mX[i,]*(vym[i,,]*res$vB3)
    ELBO_betatui[2:p]  <- ELBO_betatui[2]  + 0.5*alpha*exp(tui)*c(t(mX[i,,])%*%(vym[i,]*res$vB3))
    #apply(midtemp,2,sum)
    
    ELBO_tumuui <- muui*exp(-tu)
    ELBO_tutui  <- 0.5*exp(tui)*exp(-tu)
    
    mH.vitheta  <- matrix(c(ELBO_betamuui,ELBO_tumuui,ELBO_tvmuui,ELBO_betatui,ELBO_tutui,ELBO_tvtui),p+2,2)
    luthetaxi[[i]] <- mH.vitheta
  }
  
  lvthetaxi <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+1]
    tvj  <- vtheta[p+2*m+2*j+2]
    
    ELBO_betamuvj <- rep(0,p)
    ELBO_tumuvj <- 0
    ELBO_tvmuvj <- 0
    
    ELBO_betatvj <- rep(0,p)
    ELBO_tutvj <- 0
    ELBO_tvtvj <- 0
    
    muui <- vtheta[p + 2*(1:m) + 1]
    tui  <- vtheta[p + 2*(1:m) + 2]
    
    vmuuv <- vyx[,j]
    res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    ELBO_betamuvj[1] <- ELBO_betamuvj[1] - alpha*sum(vym[,j]*res$vB2)
    #midtemp <- mX[,j,]*(vym[,j]*res$vB2)
    ELBO_betamuvj[2:p] <- ELBO_betamuvj[2:p] - alpha*c(t(mX[,j,])%*%(vym[,j]*res$vB2))
    #apply(midtemp,2,sum)
    #sum(mX[,j]*(vym[,j]*res$vB2))
    
    ELBO_betatvj[1]  <- ELBO_betatvj[1]  + alpha*0.5*exp(tvj)*sum(vym[,j]*res$vB3)
    ELBO_betatvj[2:p]  <- ELBO_betatvj[2:p]  + alpha*0.5*exp(tvj)*c(t(mX[,j,])%*%(vym[,j]*res$vB3))
    #apply(midtemp,2,sum)
    #sum(mX[,j]*(vym[,j]*res$vB3))
    
    ELBO_tvmuvj <- muvj*exp(-tv)
    ELBO_tvtvj  <- 0.5*exp(tvj)*exp(-tv)
    
    mH.wjtheta  <- matrix(c(ELBO_betamuvj,ELBO_tumuvj,ELBO_tvmuvj,ELBO_betatvj,ELBO_tutvj,ELBO_tvtvj),p+2,2)
    
    
    #for(i in 1:m){
    #  muui <- vtheta[p+(i-1)*2+1]
    #  tui  <- vtheta[p+(i-1)*2+2]
    
    #  vmuuv <- vbeta[1] + vbeta[2]*mX[i,j]
    #  res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
    
    #  ELBO_betamuvj <- ELBO_betamuvj - c(1,mX[i,j])*alpha*vym[i,j]*as.numeric(res$vB2)
    #  ELBO_betatvj  <- ELBO_betatvj  + c(1,mX[i,j])*alpha*vym[i,j]*0.5*exp(tvj)*as.numeric(res$vB3)
    #}
    #ELBO_tvmuvj <- muvj*exp(-tv)
    #ELBO_tvtvj  <- 0.5*exp(tvj)*exp(-tv)
    
    mH.wjtheta  <- matrix(c(ELBO_betamuvj,ELBO_tumuvj,ELBO_tvmuvj,ELBO_betatvj,ELBO_tutvj,ELBO_tvtvj),p+2,2)
    
    lvthetaxi[[j]] <- mH.wjtheta
  }
  
  Hthetaxi <- matrix(0,p+2,2*m + 2*n)
  Hthetau <- matrix(0,p+2,2*m)
  Hthetav <- matrix(0,p+2,2*n)
  
  for(i in 1:m){
    Hthetau[,(2*i-1):(2*i)] <- luthetaxi[[i]] 
  }
  for(j in 1:n){
    Hthetav[,(2*j-1):(2*j)] <- lvthetaxi[[j]] 
  }
  Hthetaxi[,1:(2*m)] <- Hthetau
  Hthetaxi[,(1:(2*n))+2*m] <- Hthetav
  
  res.mH <- list(mH.00=mH.theta,lmH.0i=Hthetaxi,lmH.ii=mH.xiuv)
  
  return(list(res.vg=res.vg,res.mH=res.mH))
}


GVA.UiVj.NRDIR <- function(alpha,vtheta,family,vy,mX,gh) 
{
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3] + 1 
  
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  ###########################################################################
  ###########################################################################
  
  # Calculate 1d integrals and some auxiliary information
  res <- GVA.UiVj.vgANDmH(alpha,vtheta,family,vy,mX,gh) 
  res.vg <- res$res.vg
  res.mH <- res$res.mH
  vg.til <- res.vg$vg.0
  suv <- res.vg$vg[-c(1:(p+2))]
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
GVA.UiVj.FIT <- function(alpha,vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id) 
{
  MAXITER  <- 20
  MAXJTER  <- 20
  EPS.TERM.PAR  <- 1.0E-5
  EPS.TERM.GRAD <- 1.0E-9
  
  if (family=="LOGISTIC") {
    # Generate Quadrature Points
    N <- 3
    gh.temp <- gauss.quad(N,kind="hermite")
    #gh <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2)) 
    gh1 <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2)) 
    myGrid <- createNIGrid(dim=2, type="GHe", level=N)
    ww <- matrix(0,N^2,2)
    for(i in 1:N){
      for(j in 1:N){
        t <- (i-1)*N+j
        ww[t,2] <- sqrt(2)*gh1$w.til[i]
        ww[t,1] <- sqrt(2)*gh1$w.til[j]
      }
    }
    gh <- list(x=myGrid$nodes,w=ww)
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
  vtheta <- c(vbeta,tuv,uxi,vxi)
  sl <- GVA.UiVj.f(alpha,vtheta,family,vy,mX,gh) 
  
  cat("Initial variational log-likelihood = ",sl)
  cat("\n\n")
  
  cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
  
  for (ITER in 1:MAXITER) {
    vtheta.old <- vtheta
    
    # Calculate Newton-Raphson Direction 
    res <- GVA.UiVj.NRDIR(alpha,vtheta,family,vy,mX,gh)    
    
    stepsize <- 1                
    vtheta.temp <- vtheta + stepsize*res$dvtheta
    sl.temp <- GVA.UiVj.f(alpha,vtheta.temp,family,vy,mX,gh) 
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
        vsl[JTER] <- GVA.UiVj.f(alpha,vtheta.temp,family,vy,mX,gh) 
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
    
    if (exp(vtheta[p+1])<5.0E-9) { break; } # Break if sigma2 is getting too small
    if (exp(vtheta[p+2])<5.0E-9) { break; }
  }
  
  vbeta <- vtheta[1:p]
  sigmauv <- vtheta[p+1:2]
  
  vmuu   <- vtheta[p+2*(1:m)+1]
  vzetau <- vtheta[p+2*(1:m)+2]
  
  vmuv<- vtheta[p + 2*m + 2*(1:n) + 1]
  vzetav <- vtheta[p + 2*m + 2*(1:n) + 2]
  
  return(list(vbeta=vbeta,sigma2=exp(sigmauv),muu=vmuu,lambdau=exp(vzetau),muv=vmuv,lambdav=exp(vzetav),vbeta.serr=sqrt(diag(res$sI[1:(p+2),1:(p+2)])),sl=sl))
}


###############################################################################

###############################################################################
