###############################################################################

b0 <- function(x) { log(1+exp(x)); }
b1 <- function(x) { 1/(1+exp(-x)); }
b2 <- function(x) { mu <- b1(x); mu*(1-mu); }
b3 <- function(x) { mu <- b1(x); (mu-2*(mu^2))*(1 - mu); }
b4 <- function(x) { mu <- b1(x); (mu - 6*(mu^2)  + 6*(mu^3))*(1 - mu); }
b5 <- function(x) { mu <- b1(x); (mu - 14*(mu^2) + 36*(mu^3)  - 24*(mu^4))*(1 - mu); }
b6 <- function(x) { mu <- b1(x); (mu - 30*(mu^2) + 150*(mu^3) - 240*(mu^4) + 120*(mu^5))*(1 - mu); }

##############################################################################

B0.fun <- function(family,vmuuv,muui,tui,muvj,tvj,gh)
{ 
  if (family=="NORMALTRY"){
    vB0 <- 0.5*((vmuuv + muui)^{2} + (vmuuv + muvj)^{2} + exp(tui) + exp(tvj))
  }
  if (family=="NORMAL") {
    vB0 <- 0.5*(vmuuv + muui + muvj)^{2} + 0.5*exp(tui) +0.5*exp(tvj) 
  }
  if (family=="POISSON") {
    vB0 <- exp(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
  }
  if (family=="LOGISTIC") {
    vB0 <- log(1+exp(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj)))
    #vB0 <- B0.aghq(vmuuv,muui,tui,muvj,tvj,gh)$vB0
  }
  return(vB0)
}

###############################################################################

B12.fun <- function(family,vmuuv,muui,tui,muvj,tvj,gh) 
{
  if (family=="NORMALTRY"){
    vB1 <- 2*vmuuv + muui + muvj
    vB2 <- 1
  }
  if (family=="NORMAL"){
    vB1 <- vmmuv + muui + muvj
    vB2 <- 1
  }
  if (family=="POISSON") {
    vB1 <- exp(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
    vB2 <- vB1
  }
  if (family=="LOGISTIC") {
    vB1 <- b1(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
    vB2 <- b2(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))	
    #vB3 <- b3(mu + 0.5*sigma2)
    #vB4 <- b4(mu + 0.5*sigma2)				
    #res <- B12.aghq(vmuuv,muui,tui,muvj,tvj,gh)
    #vB1 <- res$vB1
    #vB2 <- res$vB2	
  }		
  return(list(vB1=vB1,vB2=vB2))
}

###############################################################################

B1234.fun <- function(family,vmuuv,muui,tui,muvj,tvj,gh) 
{
  
  if (family=="NORMALTRY"){
    vB1 <- 2*vmuuv + muui + muvj
    vB2 <- 1
    vB3 <- 0
    vB4 <- 0
  }
  if (family=="NORMAL") {
    vB1 <- vmuuv + muui + muvj
    vB2 <- 1
    vB3 <- 0
    vB4 <- 0
  }
  if (family=="POISSON") {
    vB1 <- exp(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
    vB2 <- vB1
    vB3 <- vB1
    vB4 <- vB1                
  } 
  if (family=="LOGISTIC") {
    vB1 <- b1(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
    vB2 <- b2(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))	
    vB3 <- b3(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))
    vB4 <- b4(vmuuv + muui + 0.5*exp(tui) + muvj + 0.5*exp(tvj))				
    #log(1+exp(mu + 0.5*sigma2)) 
    #res <- B1234.aghq(vmuuv,muui,tui,muvj,tvj,gh)
    #vB1 <- res$vB1
    #vB2 <- res$vB2	
    #vB3 <- res$vB3
    #vB4 <- res$vB4					
  }		
  return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))
}

###############################################################################
###############################################################################

find.modes <- function(vmu,tui,tvj,uv.hat) 
{
  MAXITER  <- 100
  EPS.TERM <- 1.0e-5
  
  vsigmau <- sqrt(exp(tui))
  vsigmav <- sqrt(exp(tvj))
  vsigma <- as.vector(c(vsigmau,vsigmav))
  for (ITER in 1:MAXITER) {
    uv.til <- vmu + sum(vsigma*uv.hat)
    vb0 <- b0(uv.til)
    vb1 <- b1(uv.til)
    vb2 <- b2(uv.til)
    
    vg1  <- vsigma*(vb1/vb0)  - uv.hat
    mHuu  <- (vsigmau^2)*(vb2/vb0 - (vb1/vb0)^2) - 1
    mHvv  <- (vsigmav^2)*(vb2/vb0 - (vb1/vb0)^2) - 1
    mHuv  <- vsigmau*vsigmav*(vb2/vb0 - (vb1/vb0)^2)
    
    mH <- matrix(c(mHuu,mHuv,mHuv,mHvv),2,2)
    dx  <- vg1%*%solve(mH)
    uv.hat <- uv.hat - dx
    if (all(abs(vg1)<EPS.TERM)){ break;}
  }		
  return(list(uv.hat=uv.hat,mH=mH))
}

###############################################################################

aghq.grid <- function(vmu,tui,tvj,gh) 
{
  N <- dim(gh$x)[1]
  
  sigmau <- sqrt(exp(tui))
  sigmav <- sqrt(exp(tvj))
  vsigma <- c(sigmau,sigmav)
  
  #u.0 <- faster.logistic.normal.approx.mode(vmu,exp(tui))      
  #u.0 <- (u.0 - vmu)/sigmau
  
  #v.0 <- faster.logistic.normal.approx.mode(vmu,exp(tvj))      
  #v.0 <- (v.0 - vmu)/sigmav
  
  #uv.hat <- c(u.0,v.0)
  
  uv.hat <- as.vector(c(0,0))
  #Polish the Initial Guess using Newton's method
  
  res <- find.modes(vmu,tui,tvj,uv.hat) 
  
  vmu.star     <- as.vector(res$uv.hat)
  vsigma2.star <- chol(solve(-res$mH))
  
  absvsigma2.star <- det(vsigma2.star)  
  
  vmu.til <- vmu + sum(vsigma*vmu.star)
  
  # Calculate the weights and abscissae
  vmX <- matrix(0,2,N)
  vmX.til <- matrix(0,1,N)
  vmW.til <- matrix(0,2,N)
  
  for(i in 1:N){
    vmX[,i] <- matrix(vmu.star,2,1) + vsigma2.star%*%gh$x[i,]
    vmX.til[,i] <- vmu.til + vsigma%*%(vsigma2.star%*%gh$x[i,])
    vmW.til[,i] <- sqrt(absvsigma2.star)*gh$w[i,]*exp(-0.5*vmX[,i]^2)/sqrt(2*pi)
  }  
  
  return(list(vmXstar=vmX,vmX=vmX.til,vmW=vmW.til))
}

###############################################################################

B0.aghq <- function(vmuuv,muui,tui,muvj,tvj,gh) 
{ 
  vmu <- vmuuv + muui + muvj 
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  # Calculate the weights and abscissae
  res <- aghq.grid(vmu,tui,tvj,gh) 
  
  # Evalute the integrals
  vB0 <- (res$vmW[1,]*res$vmW[2,]*b0(res$vmX))%*%matrix(1,dim(gh$w)[1],1)
  
  vB0[negind] <- vB0[negind] - vmu[negind]
  
  return(list(vB0=vB0))    
}

###############################################################################

B12.aghq <- function(vmuuv,muui,tui,muvj,tvj,gh) 
{ 
  vmu <- vmuuv + muui + muvj 
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  # Calculate the weights and abscissae
  res <- aghq.grid(vmu,tui,tvj,gh) 	
  
  # Evalute the integrals
  vB1 <- (res$vmW[1,]*res$vmW[2,]*b1(res$vmX))%*%matrix(1,dim(gh$w)[1],1)
  vB2 <- (res$vmW[1,]*res$vmW[2,]*b2(res$vmX))%*%matrix(1,dim(gh$w)[1],1)
  
  vB1[negind] <- 1 - vB1[negind]
  
  return(list(vB1=vB1,vB2=vB2))    
}

###############################################################################

B1234.aghq <- function(vmuuv,muui,tui,muvj,tvj,gh) 
{
  vmu <- vmuuv + muui + muvj 
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  res <- aghq.grid(vmu,tui,tvj,gh) 	
  
  # Evalute the integrals
  vB1 <- (res$vmW[1,]*res$vmW[2,]*b1(res$vmX))%*%matrix(1,dim(gh$w)[1],1)
  vB2 <- (res$vmW[1,]*res$vmW[2,]*b2(res$vmX))%*%matrix(1,dim(gh$w)[1],1)	
  vB3 <- (res$vmW[1,]*res$vmW[2,]*b3(res$vmX))%*%matrix(1,dim(gh$w)[1],1)	
  vB4 <- (res$vmW[1,]*res$vmW[2,]*b4(res$vmX))%*%matrix(1,dim(gh$w)[1],1)	
  
  vB1[negind] <- 1 - vB1[negind]
  vB3[negind] <-   - vB3[negind]
  
  return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))    
}

###############################################################################
###############################################################################
###############################################################################

B0rc.fun <- function(family,mu,sigma2,gh) 
{ 
  if (family=="NORMAL"){
    vB0 <- 0.5*((mu^2) + sigma2)  
  }
  if (family=="POISSON") {
    vB0 <- exp(mu+0.5*sigma2)
  } 
  if (family=="LOGISTIC") {
    vB0 <- log(1+exp(mu + 0.5*sigma2))
  }
  return(vB0)
}

###############################################################################

B12rc.fun <- function(family,mu,sigma2,gh) 
{
  if (family=="NORMAL") {
    vB1 <- mu
    vB2 <- 1
  }
  if (family=="POISSON") {
    vB1 <- exp(matrix(mu)+0.5*matrix(sigma2))
    vB2 <- vB1
  } 
  if (family=="LOGISTIC") {
    res <- B12rc.aghq(mu,sigma2,gh)
    vB1 <- res$vB1
    vB2 <- res$vB2	
  }		
  return(list(vB1=vB1,vB2=vB2))
}

###############################################################################

B1234rc.fun <- function(family,mu,sigma2,gh) 
{
  if (family=="NORMAL") {
    vB1 <- mu
    vB2 <- 1
    vB3 <- 0
    vB4 <- 0
  }
  if (family=="POISSON") {
    vB1 <- exp(mu+0.5*sigma2)
    vB2 <- vB1
    vB3 <- vB1
    vB4 <- vB1                
  } 
  if (family=="LOGISTIC"){
    vB1 <- b1(mu + 0.5*sigma2)
    vB2 <- b2(mu + 0.5*sigma2)	
    vB3 <- b3(mu + 0.5*sigma2)
    vB4 <- b4(mu + 0.5*sigma2)					
  }		
  return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))
}

###############################################################################

solve.cubic <- function(a,b,c) {
  
  Q <- (a^2 - 3*b)/9
  R <- (2*(a^3) - 9*a*b + 27*c)/54
  
  if ((R^2)<=(Q^3)) {
    theta <- acos(R/sqrt(Q^3))
    x1 <- -2*sqrt(Q)*cos(theta/3) - a/3
    x2 <- -2*sqrt(Q)*cos((theta+2*pi)/3) - a/3
    x3 <- -2*sqrt(Q)*cos((theta-2*pi)/3) - a/3
    return(c(x1,x2,x3))
  } else {
    A <- - sign(R)*(abs(R) + sqrt(R^2 - Q^3))^(1/3)
    if (A!=0) {
      B <- Q/A
    } else {
      B <- 0
    }
    return( A + B - a/3 )
  }        
}

###############################################################################

logistic.normal.approx.mode <- function(vmu,vsigma2)
{    
  a1 <- 0.0128663
  p0 <- 0.353417
  p1 <- 0.918588
  p2 <- 0.0179221
  
  vc2 <-  a1*vsigma2 + p1/p2 - vmu
  vc1 <-  a1*p1*vsigma2/p2 - 2*vsigma2 - vmu*p1/p2 + p0/p2
  vc0 <- -p1*vsigma2/p2 + a1*p0*vsigma2/p2 - vmu*p0/p2
  
  vx <- solve.cubic(vc2,vc1,vc0)
  
  return(vx)
}

###############################################################################

faster.logistic.normal.approx.mode <- function(vmu,vsigma2)
{    
  a1 <- 0.0128663
  p0 <- 0.353417
  p1 <- 0.918588
  p2 <- 0.0179221
  
  vc2 <-  a1*vsigma2 + p1/p2 - vmu
  vc1 <-  a1*p1*vsigma2/p2 - 2*vsigma2 - vmu*p1/p2 + p0/p2
  vc0 <- -p1*vsigma2/p2 + a1*p0*vsigma2/p2 - vmu*p0/p2
  
  # Almost always the middle root of the cubic solver.    
  vQ <- (vc2^2 - 3*vc1)/9
  vR <- (2*(vc2^3) - 9*vc2*vc1 + 27*vc0)/54    
  vtheta <- acos(vR/sqrt(vQ^3))
  vx <- -2*sqrt(vQ)*cos((vtheta+2*pi)/3) - vc2/3
  
  return(vx)
}

###############################################################################

findrc.modes <- function(vmu,vsigma2,vx.hat) 
{
  MAXITER  <- 100
  EPS.TERM <- 1.0e-5
  vsigma <- sqrt(vsigma2)
  for (ITER in 1:MAXITER) {
    vx.til <- vmu + vsigma*vx.hat
    vb0 <- b0(vx.til)
    vb1 <- b1(vx.til)
    vb2 <- b2(vx.til)		   
    vg  <- vsigma*(vb1/vb0)  - vx.hat
    mH  <- vsigma2*(vb2/vb0 - (vb1/vb0)^2) - 1
    dx  <- vg/mH				   
    vx.hat <- vx.hat - dx
    if (all(abs(vg)<EPS.TERM)) { break; }
  }		
  return(list(vx.hat=vx.hat,mH=mH))
}

###############################################################################

aghqrc.grid <- function(vmu,vsigma2,gh) 
{
  n <- length(vmu)
  N <- length(gh$x)
  
  vsigma <- sqrt(vsigma2)
  
  # Calculate Initial Guess  
  #vx.0 <- matrix(0,n,1)        
  #for (i in 1:n) {    
  #   # Calculate intitial guesses for modes
  #   vx.hats <- logistic.normal.approx.mode(vmu[i],vsigma2[i])  
  #   vf.hats <- b0(vx.hats)*dnorm(vx.hats,vmu[i],vsigma2[i])  
  #   # Select the intitial guess with the highest value
  #   print(c(vmu[i],vsigma2[i],vf.hats))
  #   vx.0[i] <- vx.hats[which.max(vf.hats)]
  #}
  
  
  vx.0 <- faster.logistic.normal.approx.mode(vmu,vsigma2)      
  vx.0 <- (vx.0 - vmu)/vsigma
  
  # Polish the Initial Guess using Newton's method
  res <- findrc.modes(vmu,vsigma2,vx.0) 
  vmu.star     <- as.vector(res$vx.hat)
  vsigma2.star <- as.vector(-1/res$mH)
  
  sqrt.2vsigma2.star <- as.vector(sqrt(2*vsigma2.star))
  vmu.til <- vmu + vsigma*vmu.star
  
  # Calculate the weights and abscissae
  mX     <- matrix(vmu.star,n,N) + sqrt.2vsigma2.star%o%gh$x
  mX.til <- matrix(vmu.til,n,N)  + (sqrt.2vsigma2.star*vsigma)%o%gh$x
  mW.til <- (sqrt.2vsigma2.star%o%(gh$w.til/sqrt(2*pi)))*exp(-0.5*mX^2)		
  
  return(list(mX=mX.til,mW=mW.til))
}

###############################################################################

B0rc.aghq <- function(vmu,vsigma2,gh) 
{
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  # Calculate the weights and abscissae
  res <- aghqrc.grid(as.vector(vmu),as.vector(vsigma2),gh) 
  
  # Evalute the integrals
  vB0 <- (res$mW*b0(res$mX))%*%matrix(1,length(gh$x),1)
  
  # Fix up values corresponding to originally negative vmu
  vB0[negind] <- vB0[negind] - vmu[negind]
  
  return(list(vB0=vB0))    
}

###############################################################################

B12rc.aghq <- function(vmu,vsigma2,gh) 
{
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  # Calculate the weights and abscissae
  res <- aghqrc.grid(as.vector(vmu),as.vector(vsigma2),gh) 	
  
  # Evalute the integrals
  vB1 <- (res$mW*b1(res$mX))%*%matrix(1,length(gh$x),1)
  vB2 <- (res$mW*b2(res$mX))%*%matrix(1,length(gh$x),1)
  
  # Fix up values corresponding to originally negative vmu
  vB1[negind] <- 1 - vB1[negind]	
  
  return(list(vB1=vB1,vB2=vB2))    
}

###############################################################################

B1234rc.aghq <- function(vmu,vsigma2,gh) 
{
  negind <- (vmu<0)
  vmu    <- abs(vmu)
  
  # Calculate the weights and abscissae
  res <- aghqrc.grid(as.vector(vmu),as.vector(vsigma2),gh) 	
  
  # Evalute the integrals
  vB1 <- (res$mW*b1(res$mX))%*%matrix(1,length(gh$x),1)
  vB2 <- (res$mW*b2(res$mX))%*%matrix(1,length(gh$x),1)	
  vB3 <- (res$mW*b3(res$mX))%*%matrix(1,length(gh$x),1)	
  vB4 <- (res$mW*b4(res$mX))%*%matrix(1,length(gh$x),1)	
  
  # Fix up values corresponding to originally negative vmu
  vB1[negind] <- 1 - vB1[negind]	
  vB3[negind] <-   - vB3[negind]	
  
  return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))    
}

###############################################################################

