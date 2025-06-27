
#library(statmod)
#source("CalculateB.Rs")

###############################################################################

GVA.UiVj.f <- function(vtheta,family,vy,mX,gh) 
{
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  
  #vyx0 <- matrix(0,m,n)
  #for(i in 1:m){
  #  for(j in 1:n){
  #    vyx0[i,j] <- sum(c(1,mX[i,j,])*vbeta)
  #  }
  #}
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  vyx0 <- matrix(vbeta[1],m,n) + vyx0
  
  
  elbo <- 0
  elbo_1 <- 0
  elbo_2 <- 0
  elbo_3 <- 0
  elbo_4 <- 0
  elbo_5 <- 0
  
  for (i in 1:m) {
    muui <- vtheta[p + 2*i + 1]
    tui  <- vtheta[p + 2*i + 2]
    for(j in 1:n){
      muvj <- vtheta[p + 2*m + 2*j + 1]
      tvj  <- vtheta[p + 2*m + 2*j + 2]
      
      vmuuv <- vyx0[i,j]
      vB0 <- B0.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      elbo_1 <- elbo_1 + vym[i,j]*(vmuuv + muui + muvj) - vB0
      
    }
    
    elbo_2 <- elbo_2 - 0.5*log(2*pi) - 0.5*tu - 0.5*(muui^2 + exp(tui))*exp(-tu)
    elbo_3 <- elbo_3 + 0.5*log(2*pi) + 0.5*tui + 0.5
  }
  
  for(j in 1:n){
    muvj <- vtheta[p + 2*m + 2*j + 1]
    tvj  <- vtheta[p + 2*m + 2*j + 2]
    
    elbo_4 <- elbo_4 - 0.5*log(2*pi) - 0.5*tv - 0.5*(muvj^2 + exp(tvj))*exp(-tv)
    elbo_5 <- elbo_5 + 0.5*log(2*pi) + 0.5*tvj + 0.5
  }
  
  elbo <- elbo_1 + elbo_2 + elbo_3 + elbo_4 + elbo_5
  
  if (family=="POISSON"){
    for(i in 1:m){
      for(j in 1:n){
        elbo <- elbo - lgamma(vym[i,j]+1)
      }
    }
  }
  return(elbo)
}

###############################################################################

GVA.UiVj.vg <- function(vtheta,family,vy,mX,gh) 
{ 
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  
  #vyx0 <- matrix(0,m,n)
  #for(i in 1:m){
  #  for(j in 1:n){
  #    vyx0[i,j] <- sum(c(1,mX[i,j,])*vbeta)
  #  }
  #}
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  vyx0 <- matrix(vbeta[1],m,n) + vyx0
  
  ###########################################################################
  ELBO_beta <- rep(0,p)
  ELBO_tu <- 0
  ELBO_tv <- 0
  
  #for each i
  lvg.uxi <- list()
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    ELBO_muui <- 0
    ELBO_tui  <- 0
    
    for(j in 1:n){
      muvj <- vtheta[p+2*m+2*j+1]
      tvj  <- vtheta[p+2*m+2*j+2]
      vmuuv <- vyx0[i,j]
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      ELBO_beta <- ELBO_beta + as.numeric(vym[i,j]-res$vB1)*c(1,mX[i,j,])
      ELBO_muui <- ELBO_muui + vym[i,j] - res$vB1
      ELBO_tui <- ELBO_tui - 0.5*exp(tui)*res$vB2
    }
    ELBO_tu   <- ELBO_tu   - 0.5 + 0.5*(muui^2 + exp(tui))*exp(-tu)
    ELBO_muui <- ELBO_muui - muui*exp(-tu)
    ELBO_tui  <- ELBO_tui  - 0.5*exp(tui)*exp(-tu) + 0.5
    
    lvg.uxi[[i]] <- c(ELBO_muui,ELBO_tui)
  }
  
  lvg.vxj <- list()
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+1]
    tvj  <- vtheta[p+2*m+2*j+2]
    
    ELBO_muvj <- 0
    ELBO_tvj  <- 0
    
    for (i in 1:m) {
      muui <- vtheta[p+2*i+1]
      tui  <- vtheta[p+2*i+2]
      
      vmuuv <- vyx0[i,j]
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_muvj <- ELBO_muvj + vym[i,j] - res$vB1
      ELBO_tvj <- ELBO_tvj - 0.5*exp(tvj)*res$vB2
    }
    ELBO_tv   <- ELBO_tv   - 0.5 + 0.5*(muvj^2 + exp(tvj))*exp(-tv)
    ELBO_muvj <- ELBO_muvj - muvj*exp(-tv)
    ELBO_tvj  <- ELBO_tvj  - 0.5*exp(tvj)*exp(-tv) + 0.5
    
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
GVA.UiVj.vgANDmH <- function(vtheta,family,vy,mX,gh){ 
  
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  
  #vyx0 <- matrix(0,m,n)
  #for(i in 1:m){
  #  for(j in 1:n){
  #    vyx0[i,j] <- sum(c(1,mX[i,j,])*vbeta)
  #  }
  #}
  vyx0 <- matrix(0,m,n)
  for(i in 1:(p-1)){
    vyx0 <- vyx0 + mX[,,i]*vbeta[i+1]
  }
  vyx0 <- matrix(vbeta[1],m,n) + vyx0
  
  ###########################################################################
  
  res.vg <- GVA.UiVj.vg(vtheta,family,vy,mX,gh)
  
  ###########################################################################
  mH.theta <- matrix(0,p+2,p+2)
  ELBO_betabeta <- matrix(0,p,p)
  ELBO_tutu <- 0
  ELBO_tvtv <- 0
  
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    for(j in 1:n){
      muvj <- vtheta[p+2*m+2*j+1]
      tvj  <- vtheta[p+2*m+2*j+2]
      vmuuv <- vyx0[i,j]
      
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_betabeta <- ELBO_betabeta - c(1,mX[i,j,])%*%t(c(1,mX[i,j,]))*(as.numeric(res$vB2))
    }
    ELBO_tutu <- ELBO_tutu - 0.5*(muui^2 + exp(tui))*exp(-tu)
  }
  
  for(j in 1:n){
    muvj <- vtheta[p+2*m+2*j+1]
    tvj  <- vtheta[p+2*m+2*j+2]
    
    ELBO_tvtv <- ELBO_tvtv - 0.5*(muvj^2 + exp(tvj))*exp(-tv)
  }
  
  mH.theta[1:p,1:p] <- ELBO_betabeta 
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
    
    for(j in 1:n){
      muvj <- vtheta[p+2*m+2*j+1]
      tvj  <- vtheta[p+2*m+2*j+2]
      
      vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_muuimuui <- ELBO_muuimuui - (res$vB2)
      ELBO_muuitui  <- ELBO_muuitui  - (res$vB3)*exp(tui)*0.5
      ELBO_tuitui   <- ELBO_tuitui   - (res$vB4)*exp(2*tui)*0.25 - (res$vB2)*exp(tui)*0.5 
    }
    
    ELBO_muuimuui <- ELBO_muuimuui - exp(-tu) 
    ELBO_tuitui   <- ELBO_tuitui   - exp(-tu)*exp(tui)*0.5 
    
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
    
    for(i in 1:m){
      muui <- vtheta[p+2*i+1]
      tui  <- vtheta[p+2*i+2]
      vmuuv <- vyx0[i,j]
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_muvjmuvj <- ELBO_muvjmuvj - (res$vB2)
      ELBO_muvjtvj  <- ELBO_muvjtvj  - (res$vB3)*exp(tvj)*0.5
      ELBO_tvjtvj   <- ELBO_tvjtvj   - (res$vB4)*exp(2*tvj)*0.25 - (res$vB2)*exp(tvj)*0.5
    }
    ELBO_muvjmuvj <- ELBO_muvjmuvj - exp(-tv) 
    ELBO_tvjtvj   <- ELBO_tvjtvj   - 0.5*exp(tvj)*exp(-tv) 
    
    mH.vjvj <- matrix(c(ELBO_muvjmuvj,ELBO_muvjtvj,ELBO_muvjtvj,ELBO_tvjtvj),2,2) 
    
    lwxi.vv[[j]] <- mH.vjvj
    
  }
  
  lwxi.uivj <- list()
  for(i in 1:m){
    muui <- vtheta[p+2*i+1]
    tui  <- vtheta[p+2*i+2]
    
    for(j in 1:n){
      muvj <- vtheta[p+2*m+2*j+1]
      tvj  <- vtheta[p+2*m+2*j+2]
      vmuuv <- vyx0[i,j]
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_muuimuvj <- -(res$vB2)
      ELBO_muuitvj  <- -(res$vB3)*0.5*exp(tvj)
      ELBO_tuimuvj  <- -(res$vB3)*0.5*exp(tui)
      ELBO_tuitvj   <- -(res$vB4)*0.25*exp(tui+tvj)
      
      index <- (i-1)*n + j
      
      lwxi.uivj[[index]] <- matrix(c(ELBO_muuimuvj,ELBO_tuimuvj,ELBO_muuitvj,ELBO_tuitvj),2,2)
      
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
    
    for(j in 1:n){
      muvj <- vtheta[p+2*m+2*j+1]
      tvj  <- vtheta[p+2*m+2*j+2]
      vmuuv <- vyx0[i,j]
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_betamuui <- ELBO_betamuui - c(1,mX[i,j,])*as.numeric(res$vB2)
      ELBO_betatui  <- ELBO_betatui  - c(1,mX[i,j,])*0.5*exp(tui)*as.numeric(res$vB3)
    }
    
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
    
    for(i in 1:m){
      muui <- vtheta[p+2*i+1]
      tui  <- vtheta[p+2*i+2]
      #vmuuv <- vbeta[1] + sum(vbeta[2:p]*mX[i,j,])
      vmuuv <- vyx0[i,j]
      res <- B1234.fun(family,vmuuv,muui,tui,muvj,tvj,gh)
      
      ELBO_betamuvj <- ELBO_betamuvj - c(1,mX[i,j,])*as.numeric(res$vB2)
      ELBO_betatvj  <- ELBO_betatvj  - c(1,mX[i,j,])*0.5*exp(tvj)*as.numeric(res$vB3)
    }
    ELBO_tvmuvj <- muvj*exp(-tv)
    ELBO_tvtvj  <- 0.5*exp(tvj)*exp(-tv)
    
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


GVA.UiVj.NRDIR <- function(vtheta,family,vy,mX,gh) 
{
  m <- dim(mX)[1]
  n <- dim(mX)[2]
  p <- dim(mX)[3]+1
  ###########################################################################
  vbeta <- vtheta[1:p]
  tu <- vtheta[p+1]
  tv <- vtheta[p+2]
  vym <- vy
  ###########################################################################
  
  # Calculate 1d integrals and some auxiliary information
  res <- GVA.UiVj.vgANDmH(vtheta,family,vy,mX,gh) 
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
GVA.UiVj.FIT <- function(vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id) 
{
  MAXITER  <- 100
  MAXJTER  <- 20
  EPS.TERM.PAR  <- 5.0E-4
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
  vbeta <- rep(0,p)
  vtheta <- c(vbeta,tuv,uxi,vxi)
  
  sl <- GVA.UiVj.f(vtheta,family,vy,mX,gh) 
  NI <- 0
  
  cat("Initial variational log-likelihood = ",sl)
  cat("\n\n")
  
  cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
  for (ITER in 1:MAXITER) {
    #ITER <- 2
    vtheta.old <- vtheta
    
    # Calculate Newton-Raphson Direction 
    res <- GVA.UiVj.NRDIR(vtheta,family,vy,mX,gh)    
    
    stepsize <- 1                
    vtheta.temp <- vtheta + stepsize*res$dvtheta
    sl.temp <- GVA.UiVj.f(vtheta.temp,family,vy,mX,gh) 
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
        vsl[JTER] <- GVA.UiVj.f(vtheta.temp,family,vy,mX,gh) 
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
    if (err.vg<EPS.TERM.GRAD) { break; }
    if (err.par<EPS.TERM.PAR) { break; }
    
    if (exp(vtheta[p+1])<EPS.TERM.PAR) { break; } # Break if sigma2 is getting too small
    if (exp(vtheta[p+2])<EPS.TERM.PAR) { break; }
  }
  
  vbeta   <- vtheta[1:p]
  sigmauv <- vtheta[p+1:2]
  
  
  vmuu   <- vtheta[p+1+2*(1:m)]
  vzetau <- vtheta[p+2+2*(1:m)]
  
  vmuv   <- vtheta[p + 1 + 2*m + 2*(1:n)]
  vzetav <- vtheta[p + 2 + 2*m + 2*(1:n)]
  
  return(list(NI = NI,vbeta=vbeta,sigma2=exp(sigmauv),muu=vmuu,lambdau=exp(vzetau),muv=vmuv,lambdav=exp(vzetav),vbeta.serr=sqrt(abs(diag(res$sI))),sl=sl))
}


###############################################################################
