########## Simulation for Logistic case ##########
# Last changed: 20 Aug 2024

###############################################################################
rm(list = ls())
setwd("C:\\Users\\DELL\\Desktop\\Resubmit-SS-R3\\Logistic")
#setwd("C:\\Users\\DELL\\Desktop\\Resubmit-SS-R2\\R code for GVACL\\Logistic")
library(Matrix)
library(MASS)
library(statmod)
library(glmmTMB)
library(mvQuad)
library(statmod)
library(mvtnorm)
source("CalculateB.Rs")
source("CLGVA.R")
source("GVAuivj.R")
do.GVAUiVj <- TRUE
do.GVACRCL <- TRUE
do.Laplace <- TRUE
###############################################################################
###Simulation settings###
Sigma <- c(0.5,0.5)

delta <- 1.1
M <- c(200,250,300,400)
N <- round(M^(1/delta),0)
beta.true <-  c(-2,1,-0.5,1,-1,0.5,0.8,-0.4,0.1)
p <- length(beta.true)

TRIALS  <- 500  # Total number of trials
family <- c("LOGISTIC")
###############################################################################
### Record results###
theta.gva          <- matrix(0,length(M)*TRIALS*(p+2),1); dim(theta.gva)          <- c(length(M),TRIALS,p+2)
theta.regva        <- matrix(0,length(M)*TRIALS*(p+3),1); dim(theta.regva)        <- c(length(M),TRIALS,p+3)
theta.lap          <- matrix(0,length(M)*TRIALS*(p+2),1); dim(theta.lap)          <- c(length(M),TRIALS,p+2)

se.theta.gva       <- matrix(0,length(M)*TRIALS*(p+2),1); dim(se.theta.gva)       <- c(length(M),TRIALS,p+2)
se.theta.regva     <- matrix(0,length(M)*TRIALS*(p+3),1); dim(se.theta.regva)     <- c(length(M),TRIALS,p+3)
asy.se.theta.regva <- matrix(0,length(M)*TRIALS*(p+2),1); dim(asy.se.theta.regva) <- c(length(M),TRIALS,p+2)
se.theta.lap       <- matrix(0,length(M)*TRIALS*(p+2),1); dim(se.theta.lap)       <- c(length(M),TRIALS,p+2)

times.gva    <- matrix(0,length(M)*TRIALS,1); dim(times.gva) <- c(length(M),TRIALS)
times.regva  <- matrix(0,length(M)*TRIALS,1); dim(times.regva) <- c(length(M),TRIALS)
times.lap    <- matrix(0,length(M)*TRIALS,1); dim(times.lap) <- c(length(M),TRIALS)

converged.gva    <- matrix(0,length(M)*TRIALS); dim(converged.gva) <- c(length(M),TRIALS)
converged.regva  <- matrix(0,length(M)*TRIALS); dim(converged.regva) <- c(length(M),TRIALS)
converged.lap    <- matrix(0,length(M)*TRIALS); dim(converged.lap) <- c(length(M),TRIALS)

ITER.gva <-  matrix(0,length(M)*TRIALS); dim(ITER.gva) <- c(length(M),TRIALS)
ITER.regva <-  matrix(0,length(M)*TRIALS); dim(ITER.regva) <- c(length(M),TRIALS)
ITER.lap <-  matrix(0,length(M)*TRIALS); dim(ITER.lap) <- c(length(M),TRIALS)


################################################
ITER0 <- 1
trial0 <- 1
for (ITER in ITER0:length(M)) {
  for (trial in trial0:TRIALS) {
    #############################################################################
    m  <- M[ITER]
    n  <- N[ITER]
    sigmau.true <- Sigma[1]
    sigmav.true <- Sigma[2]
    
    set.seed(trial)
    
    vuu.true <- rnorm(m,0,sigmau.true)
    vuv.true <- rnorm(n,0,sigmav.true)
    
    veta <- c()
    vx   <- c()
    id   <- c()
    Xmu <- rep(0,p-1)
    XSig <- matrix(0,p-1,p-1)
    diag(XSig) <- rep(1,p-1)
    X <- matrix(0,m*n*(p-1),1); dim(X) <- c(m,n,p-1)
    
    for (i in 1:m) {
      for (j in 1:n) {
        id <- c(id,i)
        x <- rmvnorm(1,Xmu,XSig)
        X[i,j,] <- x
        veta <- c(veta,sum(beta.true*c(1,x)) + vuu.true[i] + vuv.true[j])
      }
    }
    idj  <- c()
    for (j in 1:n) {
      for (i in 1:m) {
        idj <- c(idj,j)
      }
    }
    
    vmu.true <- 1/(1+exp(-veta))
    rp <- rbinom(n*m,1,vmu.true)
    vy <- t(matrix(rp,n,m))
    
    mX <- X
    #############################################################################
    vbeta  <- rep(0,p+1)
    vSigma <- c(1,1)
    lvmu <- list()
    lvLambda <- list()
    for(j in 1:n){
      lvmu[[j]]     <- matrix(0,1,1)
      lvLambda[[j]] <- matrix(1,1,1)
    }
    
    lumu <- list()
    luLambda <- list()
    for (i in 1:m) {
      lumu[[i]]     <- matrix(0,1,1)
      luLambda[[i]] <- matrix(1,1,1) 
    }
    
    family = "LOGISTIC"
    
    if (do.GVACRCL){
      cat("trial=",trial,"    m=",m," n=",n,"    sigmaw=",Sigma[1]," sigmav=",Sigma[2],"\n")
      cat("Fitting model using ReGVA.CLCR method \n")
      
      bval1  <- proc.time()
      res <- try(ReGVA.CLCR.fit <- ReGVA.CLCR.FIT(vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id),silent=TRUE)
      eval1  <- proc.time()
      
      if ((class(res)[1])!="try-error") {								
        
        theta.regva[ITER,trial,1:(p+1)] <- res$vbeta
        theta.regva[ITER,trial,p+2:3] <- sqrt(res$sigma2)
        times.regva[ITER,trial]  <- eval1[3]-bval1[3]	
        se.theta.regva[ITER,trial,]     <- res$vbeta.serr
        ITER.regva[ITER,trial] <- res$NI
        converged.regva[ITER,trial] <- TRUE    
      } else {
        converged.regva[ITER,trial] <- FALSE
      }
    }
    
    
    ###############################################################################
    if (do.GVAUiVj){
      source("GVAuivj.R")
      cat("trial=",trial,"    m=",m,"    sigmaw=",Sigma[1]," sigmav=",Sigma[2],"\n")
      cat("Fitting model using GVA.UiVj \n")
      
      vy <- t(matrix(rp,n,m))
      mX <- X
      #############################################################################
      vbeta  <- c(0,0)
      vSigma <- c(1,1)
      lvmu <- list()
      lvLambda <- list()
      for(j in 1:n){
        lvmu[[j]]     <- matrix(0,1,1)
        lvLambda[[j]] <- matrix(1,1,1)
      }
      
      lumu <- list()
      luLambda <- list()
      for (i in 1:m) {
        lumu[[i]]     <- matrix(0,1,1)
        luLambda[[i]] <- matrix(1,1,1) 
      }
      
      bval2  <- proc.time()
      res1 <- try(GVA.UiVj.fit <- GVA.UiVj.FIT(vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id),silent=TRUE)
      eval2  <- proc.time()
      
      if ((class(res1)[1])!="try-error") {			
        theta.gva[ITER,trial,1:p] <- res1$vbeta
        theta.gva[ITER,trial,p+(1:2)] <- sqrt(res1$sigma2)
        
        se.theta.gva[ITER,trial,1:(p+2)] <- res1$vbeta.serr
        times.gva[ITER,trial] <- eval2[3]-bval2[3]	
        ITER.gva[ITER,trial] <- res1$NI
        converged.gva[ITER,trial] <- TRUE    
      } else {
        converged.gva[ITER,trial] <- FALSE
      }	
    }
    
    ###############################################################################
    ###glmmTMB###
    if (do.Laplace){
      samp <- rep(c(1:n),m)
      lapx <- matrix(0,m*n,dim(mX)[3])
      for(i in 1:dim(mX)[3]){
        lapx[,i] <- as.vector(t(mX[,,i]))
      }
      
      dat <- data.frame(y=rp,vx=lapx,id=idj,sample=samp)
      
      bval3  <- proc.time()
      m3 <- try(Lap.fit <- glmmTMB(y ~ vx.1 + vx.2 + vx.3 + vx.4 + vx.5 + vx.6 + vx.7 + vx.8 + (1|id) + (1|sample), ziformula=~0, family=binomial, data=dat), silent=TRUE)
      eval3  <- proc.time()
      
      if ((class(m3)[1])!="try-error") {			
        
        theta.lap[ITER,trial,] <- c(m3$fit$par[1:p],sqrt(summary(m3)$varcor$cond$id[1]),sqrt(summary(m3)$varcor$cond$sample[1]))
        
        se.theta.lap[ITER,trial,1:p] <- summary(m3)$coefficients$cond[1:p,2]
        
        times.lap[ITER,trial] <- eval3[3]-bval3[3]	
        
        converged.lap[ITER,trial] <- TRUE    
      } else {
        converged.lap[ITER,trial] <- FALSE
      }	
    }
    save.image("Log-S-M500-delta-1-1.Rdata")
  }
}

###############################################################################

load("C:\\Users\\DELL\\Desktop\\Resubmit-SS-R3\\Logistic\\Log-S-M500-delta-1-1.Rdata")

indOutliers <- function(x,nsd=5) {
  z.scores <- (x - median(x))/mad(x)
  return(which(abs(z.scores)>nsd))
}

gva.beta0 <- list()
gva.beta1 <- list()
gva.sigma1 <- list()
gva.sigma2 <- list()

lap.beta0 <- list()
lap.beta1 <- list()
lap.sigma1 <- list()
lap.sigma2 <- list()

regva.beta0 <- list()
regva.beta1 <- list()
regva.sigma1 <- list()
regva.sigma2 <- list()

gva.time <- c()
lap.time <- c()
regva.time <- c()

effsize <- list()
TRMT <- 500
for (ITER in 1:length(M)) {
  m <- M[ITER]
  n <- N[ITER]
  sigma.true1 <- Sigma[1]  
  sigma.true2 <- Sigma[2]
  
  
  cat("m=",m,"  n=",n,"  beta0=",beta.true[1]," beta1=",beta.true[2]," beta2=",beta.true[3],"  sigma1=",sigma.true1,"  sigma2=",sigma.true2,"\n")
  # Remove cases where any one of the methods did not converge         
  inds <- c(1:length(converged.gva[ITER,]))
  uncong.gva   <- which(converged.gva[ITER,] !=1)
  uncong.regva <- which(converged.regva[ITER,] !=1)  
  uncong.lap   <- which(converged.lap[ITER,] !=1 )
  uncong.gva <- c(uncong.gva,which(ITER.gva[ITER,]==100))
  uncong.regva <- c(uncong.regva,which(ITER.regva[ITER,]==100))
  
  inds.outliers <- c()
  for(i in 1:dim(theta.lap)[3]){
    inds.outliers <- c(inds.outliers,
                       indOutliers(theta.regva[ITER,1:TRMT,i]),
                       indOutliers(theta.lap[ITER,1:TRMT,i]),
                       indOutliers(theta.gva[ITER,1:TRMT,i]),
                       which(asy.se.theta.regva[ITER,1:TRMT,i]== "NaN"),
                       which(se.theta.lap[ITER,1:TRMT,i]=="NaN"),
                       which(se.theta.gva[ITER,1:TRMT,i]=="NaN")
    )
  }
  inds.outliers <- c(inds.outliers,
                     indOutliers(theta.regva[ITER,1:TRMT,p+3]),
                     which(se.theta.regva[ITER,1:TRMT,p+3]=="NaN"))
  inds.outliers <- c(uncong.gva,uncong.lap,uncong.regva,inds.outliers)
  inds.outliers <- unique(inds.outliers)
  if (length(inds.outliers)>0) {              
    inds <- inds[-inds.outliers]     
  }
  
  effsize[[ITER]] <- inds
  
  print(length(inds))
  
  locfun <- mean
  sprfun <- sd
  ######################################################################
  # Calculate and print summary statistic for fitted beta0 values
  
  cat("beta0 \n")
  
  dp <- 4
  
  ### GVA ###
  mean.beta0.gva    <- round(locfun(theta.gva[ITER,inds,1]),dp)
  mean.beta0.lap    <- round(locfun(theta.lap[ITER,inds,1]),dp)
  mean.beta01.regva <- round(locfun(theta.regva[ITER,inds,1]),dp)
  mean.beta02.regva <- round(locfun(theta.regva[ITER,inds,2]),dp)
  mean.beta0.regva  <- round(locfun((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,5]^2/2-theta.regva[ITER,inds,6]^2/2)/2),dp)
  meseasy.beta0.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,1]),dp)
  
  sd.beta0.gva     <- round(sprfun(theta.gva[ITER,inds,1]),dp)  
  sd.beta0.lap     <- round(sprfun(theta.lap[ITER,inds,1]),dp)
  sd.beta0.regva   <- round(sprfun((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,5]^2/2-theta.regva[ITER,inds,6]^2/2)/2),dp)
  
  mstderr.beta0.gva    <- round(locfun(se.theta.gva[ITER,inds,1]),dp)
  mstderr.beta0.lap    <- round(locfun(se.theta.lap[ITER,inds,1]),dp)
  
  rmse.beta0.gva <-  round(sqrt(mean((theta.gva[ITER,inds,1] - beta.true[1])^{2})),dp)
  rmse.beta0.lap <-  round(sqrt(mean((theta.lap[ITER,inds,1] - beta.true[1])^{2})),dp)
  rmse.beta0.regva <- round(sqrt(locfun(((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,p+2]^2/2-theta.regva[ITER,inds,p+3]^2/2)/2 - beta.true[1])^2)),dp)
  
  cat("beta0.gva \n")
  cat(mean.beta0.gva, "(" ,sd.beta0.gva , ") & " , mstderr.beta0.gva, "&", rmse.beta0.gva, "\n") 
  cat("beta0.lap \n")
  cat(mean.beta0.lap,"(",sd.beta0.lap ,") & ", mstderr.beta0.lap,"&", rmse.beta0.lap, " \n") 
  cat("beta0.regva \n")
  cat(mean.beta0.regva, "(",sd.beta0.regva ,") & ", meseasy.beta0.regvaasy ,"&", rmse.beta0.regva, " \n")  
  
  cat("beta1 \n")
  mean.beta1.gva   <- round(locfun(theta.gva[ITER,inds,2]),dp)
  mean.beta1.lap   <- round(locfun(theta.lap[ITER,inds,2]),dp)
  mean.beta1.regva <- round(locfun(theta.regva[ITER,inds,3]),dp)
  
  sd.beta1.gva   <- round(sprfun(theta.gva[ITER,inds,2]),dp)  
  sd.beta1.lap   <- round(sprfun(theta.lap[ITER,inds,2]),dp)  
  sd.beta1.regva <- round(sprfun(theta.regva[ITER,inds,3]),dp)  
  
  mstderr.beta1.gva   <- round(locfun(se.theta.gva[ITER,inds,2]),dp)
  mstderr.beta1.lap   <- round(locfun(se.theta.lap[ITER,inds,2]),dp)
  meseasy.beta1.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,2]),dp)
  
  rmse.beta1.gva <-  round(sqrt(mean((theta.gva[ITER,inds,2] - beta.true[2])^{2})),dp)
  rmse.beta1.lap <-  round(sqrt(mean((theta.lap[ITER,inds,2] - beta.true[2])^{2})),dp)
  rmse.beta1.regva <- round(sqrt(locfun((theta.regva[ITER,inds,3] - beta.true[2])^{2})),dp)
  
  
  cat("beta1.gva \n")
  cat(mean.beta1.gva,"(",sd.beta1.gva ,") & ", mstderr.beta1.gva, "& ", rmse.beta1.gva, " \n")
  cat("beta1.lap \n")
  cat(mean.beta1.lap,"(",sd.beta1.lap ,") & ", mstderr.beta1.lap, "& ", rmse.beta1.lap, " \n")
  cat("beta1.regva \n")
  cat(mean.beta1.regva, "(",sd.beta1.regva ,") & ", meseasy.beta1.regvaasy, "& ", rmse.beta1.regva, " \n")
  
  cat("beta2 \n")
  
  mean.beta2.gva   <- round(locfun(theta.gva[ITER,inds,3]),dp)
  mean.beta2.lap   <- round(locfun(theta.lap[ITER,inds,3]),dp)
  mean.beta2.regva <- round(locfun(theta.regva[ITER,inds,4]),dp)
  
  sd.beta2.gva   <- round(sprfun(theta.gva[ITER,inds,3]),dp)  
  sd.beta2.lap   <- round(sprfun(theta.lap[ITER,inds,3]),dp)  
  sd.beta2.regva <- round(sprfun(theta.regva[ITER,inds,4]),dp)  
  
  mstderr.beta2.gva   <- round(locfun(se.theta.gva[ITER,inds,3]),dp)
  mstderr.beta2.lap  <- round(locfun(se.theta.lap[ITER,inds,3]),dp)
  meseasy.beta2.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,3]),dp)
  
  rmse.beta2.gva <-  round(sqrt(mean((theta.gva[ITER,inds,3] - beta.true[3])^{2})),dp)
  rmse.beta2.lap <-  round(sqrt(mean((theta.lap[ITER,inds,3] - beta.true[3])^{2})),dp)
  rmse.beta2.regva <- round(sqrt(mean((theta.regva[ITER,inds,4] - beta.true[3])^{2})),dp)
  
  
  cat("beta2.gva \n")
  cat(mean.beta2.gva,"(",sd.beta2.gva ,") & ", mstderr.beta2.gva, "& ", rmse.beta2.gva, " \n")
  cat("beta2.lap \n")
  cat(mean.beta2.lap,"(",sd.beta2.lap ,") & ", mstderr.beta2.lap, "& ", rmse.beta2.lap, " \n")
  cat("beta2.regva \n")
  cat(mean.beta2.regva, "(",sd.beta2.regva ,") & ", meseasy.beta2.regvaasy, "& ", rmse.beta2.regva, " \n")
  
  cat("sigma \n")
  
  mean.sigmau.gva   <- round(locfun(theta.gva[ITER,inds,p+1]),dp)
  mean.sigmav.gva   <- round(locfun(theta.gva[ITER,inds,p+2]),dp)
  
  mean.sigmau.lap   <- round(locfun(theta.lap[ITER,inds,p+1]),dp)
  mean.sigmav.lap   <- round(locfun(theta.lap[ITER,inds,p+2]),dp)
  
  mean.sigmau.regva <- round(locfun(theta.regva[ITER,inds,p+2]),dp)
  mean.sigmav.regva <- round(locfun(theta.regva[ITER,inds,p+3]),dp)
  
  sd.sigmau.gva   <- round(sprfun(theta.gva[ITER,inds,p+1]),dp)  
  sd.sigmav.gva   <- round(sprfun(theta.gva[ITER,inds,p+2]),dp)
  
  sd.sigmau.lap   <- round(sprfun(theta.lap[ITER,inds,p+1]),dp)  
  sd.sigmav.lap   <- round(sprfun(theta.lap[ITER,inds,p+2]),dp)  
  
  sd.sigmau.regva <- round(sprfun(theta.regva[ITER,inds,p+2]),dp)
  sd.sigmav.regva <- round(sprfun(theta.regva[ITER,inds,p+3]),dp)
  
  
  mstderr.sigma1.gva <-  round(locfun(se.theta.gva[ITER,inds,p+1]*theta.gva[ITER,inds,p+1]/2),dp)
  mstderr.sigma2.gva <-  round(locfun(se.theta.gva[ITER,inds,p+2]*theta.gva[ITER,inds,p+2]/2),dp)
  
  #mstderr.sigma1.lap <-  round(locfun(stderr.sigma1.lap[ITER,JTER,inds]*sigmau.lap[ITER,JTER,inds]/2),dp)
  #mstderr.sigma2.lap <-  round(locfun(stderr.sigma2.lap[ITER,JTER,inds]*sigmav.lap[ITER,JTER,inds]/2),dp)
  
  
  mstderr.sigma1.lap <-  round(locfun(se.theta.lap[ITER,inds,p+1]*theta.lap[ITER,inds,p+1]/2),dp)
  mstderr.sigma2.lap <-  round(locfun(se.theta.lap[ITER,inds,p+2]*theta.lap[ITER,inds,p+2]/2),dp)
  
  meseasy.sigma1.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,p+1]),dp)
  meseasy.sigma2.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,p+2]),dp)
  
  rmse.sigma1.gva <-  round(sqrt(mean((theta.gva[ITER,inds,p+1] - Sigma[1])^{2})),dp)
  rmse.sigma2.gva <-  round(sqrt(mean((theta.gva[ITER,inds,p+2] - Sigma[2])^{2})),dp)
  
  rmse.sigma1.lap <-  round(sqrt(mean((theta.lap[ITER,inds,p+1] - Sigma[1])^{2})),dp)
  rmse.sigma2.lap <-  round(sqrt(mean((theta.lap[ITER,inds,p+2] - Sigma[2])^{2})),dp)
  
  rmse.sigma1.regva <- round(sqrt(mean((theta.regva[ITER,inds,p+2] - Sigma[1])^{2})),dp)
  rmse.sigma2.regva <- round(sqrt(mean((theta.regva[ITER,inds,p+3] - Sigma[2])^{2})),dp)
  
  cat(mean.sigmau.gva , "(",sd.sigmau.gva ,") &  ", mstderr.sigma1.gva, " & ", rmse.sigma1.gva, " \n") 
  cat(mean.sigmav.gva , "(",sd.sigmav.gva ,") &  ", mstderr.sigma2.gva, " & ", rmse.sigma2.gva, " \n") 
  
  cat(mean.sigmau.lap , "(",sd.sigmau.lap ,") &  ", mstderr.sigma1.lap, " & ", rmse.sigma1.lap, " \n") 
  cat(mean.sigmav.lap , "(",sd.sigmav.lap ,") &  ", mstderr.sigma2.lap, " & ", rmse.sigma2.lap," \n") 
  
  cat(mean.sigmau.regva ,"(",sd.sigmau.regva ,") &  ", meseasy.sigma1.regvaasy, " & ", rmse.sigma1.regva," \n")
  cat(mean.sigmav.regva ,"(",sd.sigmav.regva ,") &  ", meseasy.sigma2.regvaasy, " & ", rmse.sigma2.regva," \n")
  
  cat("times - ", round(mean(times.gva[ITER,inds]),dp),
      round(mean(times.lap[ITER,inds]),dp),
      round(mean(times.regva[ITER,inds]),dp),"\n") 
  
  cat("\n\n") 
  
  gva.time <- c(gva.time, round(mean(times.gva[ITER,inds]),dp))
  lap.time <- c(lap.time, round(mean(times.lap[ITER,inds]),dp))
  regva.time <- c(regva.time, round(mean(times.regva[ITER,inds]),dp))
  
  ###################################################################### 
}
M
N
###Figure S3 ###
gvatime <- data.frame(x=log(M*N,10),y=log(gva.time,10))
laptime <- data.frame(x=log(M*N,10),y=log(lap.time,10))
regvatime <- data.frame(x=log(M*N,10),y=log(regva.time,10))
gvatime <- data.frame(x=M*N,y=gva.time)
laptime <- data.frame(x=M*N,y=lap.time)
regvatime <- data.frame(x=M*N,y=regva.time)
Log.laptime <- laptime 
###plot the figure of CPU Time vs sample size###
library(ggplot2)
log10_labels <- function(x) {
  parse(text = paste0("10^", log10(x)))  # 转换为 10^x 格式
}
plot1 <- 
  ggplot()+                   # 添加连线
  geom_point(data = gvatime, aes(x=x, y=y, color = "blue"), size =3) + 
  geom_point(data = laptime, aes(x=x, y=y, color = "red"), size = 3) +
  geom_point(data = regvatime, aes(x=x, y=y, color = "green"), size =3) +
  geom_line(data = gvatime, aes(x=x, y=y, color = "blue")) +
  geom_line(data = laptime, aes(x=x, y=y, color = "red")) +
  geom_line(data = regvatime, aes(x=x, y=y, color = "green")) +
  scale_color_manual(values =c("blue","green","red"),labels = c("GVA","GVACL","glmmTMB")) +
  scale_y_continuous(
    trans = "log10",                # 纵轴对数变换
    limits = c(5e-2, 5e2),          # 设置纵轴范围
    breaks = c(1e-1, 1e0, 1e1, 1e2, 1e3),  # 设置刻度
    labels = log10_labels           # 使用 10^x 格式标签
  ) +
  scale_x_continuous(
    trans = "log10",                # 横轴对数变换
    limits = c(1e2, 1e5),           # 设置横轴范围
    breaks = c(1e2, 1e3, 1e4, 1e5), # 设置刻度
    labels = log10_labels           # 使用 10^x 格式标签
  ) +
  labs(
    title = "CPU sec. Logistic",     # 图表标题
    x = "Sample size",              # 横轴标签
    y = "Time/s",                   # 纵轴标签
    color = "Method"                # 图例标题
  )+
  labs(title = "CPU sec. Logistic", x= "Sample size", y="Time/s") +
  theme_minimal() + 
  theme(legend.position = "top")

plot1 <- ggplot() + 
  geom_point(data = gvatime, aes(x=x, y=y, color = "blue"), size =3) + 
  geom_point(data = laptime, aes(x=x, y=y, color = "red"), size = 3) +
  geom_point(data = regvatime, aes(x=x, y=y, color = "green"), size =3) +
  geom_line(data = gvatime, aes(x=x, y=y, color = "blue")) +
  geom_line(data = laptime, aes(x=x, y=y, color = "red")) +
  geom_line(data = regvatime, aes(x=x, y=y, color = "green")) +
  scale_color_manual(values =c("blue","green","red"),labels = c("GVA","GVACL","glmmTMB")) +
  scale_x_continuous(breaks = log10(c(1e2,1e3,1e4,1e5)),
                     labels = scales::label_math(10^.x)) +
  scale_y_continuous(breaks = log10(c(1e-2,1e-1,1e0,1e1,1e2)),
                     labels = scales::label_math(10^.x)) +
  labs(title = "CPU sec. logistic", x= "Sample size", y="Time/s") +
  theme_minimal() + 
  theme(legend.position = "top")

data.gva0 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.gva[1,effsize[[1]],1],
            theta.gva[2,effsize[[2]],1],
            theta.gva[3,effsize[[3]],1],
            theta.gva[4,effsize[[4]],1])
)

data.gva0$group <- factor(data.gva0$group,
                          levels = c("24800","37750","53700","92800")
)
beta0gvabox <- ggplot(data.gva0, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") +
  geom_hline( yintercept = -2, color ="red") + 
  scale_y_continuous(limits = c(-2.25, -1.75)) +
  labs(title=expression(beta[1]~"GVA"), x="Sample size ", y ="" ) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.lap0 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.lap[1,effsize[[1]],1],
            theta.lap[2,effsize[[2]],1],
            theta.lap[3,effsize[[3]],1],
            theta.lap[4,effsize[[4]],1])
)

data.lap0$group <- factor(data.lap0$group,
                          levels = c("24800","37750","53700","92800")
)
beta0lapbox <- ggplot(data.lap0, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -2, color ="red") + 
  scale_y_continuous(limits = c(-2.25, -1.75)) +
  labs(title=expression(beta[1]~"glmmTMB"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.regva0 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c((theta.regva[1,effsize[[1]],1] + theta.regva[1,effsize[[1]],2] - theta.regva[1,effsize[[1]],p+2]^2/2-theta.regva[1,effsize[[1]],p+3]^2/2)/2,
            (theta.regva[2,effsize[[2]],1] + theta.regva[2,effsize[[2]],2] - theta.regva[2,effsize[[2]],p+2]^2/2-theta.regva[2,effsize[[2]],p+3]^2/2)/2,
            (theta.regva[3,effsize[[3]],1] + theta.regva[3,effsize[[3]],2] - theta.regva[3,effsize[[3]],p+2]^2/2-theta.regva[3,effsize[[3]],p+3]^2/2)/2,
            (theta.regva[4,effsize[[4]],1] + theta.regva[4,effsize[[4]],2] - theta.regva[4,effsize[[4]],p+2]^2/2-theta.regva[4,effsize[[4]],p+3]^2/2)/2
  )
)

data.regva0$group <- factor(data.regva0$group,
                            levels = c("24800","37750","53700","92800")
)
beta0regvabox <- ggplot(data.regva0, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -2, color ="red") + 
  scale_y_continuous(limits = c(-2.25, -1.75)) +
  labs(title=expression(beta[1]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))
library(patchwork)
beta0gvabox  + beta0lapbox + beta0regvabox + plot_layout(ncol = 3)

data.gva1 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.gva[1,effsize[[1]],2],
            theta.gva[2,effsize[[2]],2],
            theta.gva[3,effsize[[3]],2],
            theta.gva[4,effsize[[4]],2])
)

data.gva1$group <- factor(data.gva1$group,
                          levels = c("24800","37750","53700","92800")
)
beta1gvabox <- ggplot(data.gva1, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 1, color ="red") + 
  scale_y_continuous(limits = c(0.9, 1.1)) +
  labs(title=expression(beta[2]~"GVA"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))


data.lap1 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.lap[1,effsize[[1]],2],
            theta.lap[2,effsize[[2]],2],
            theta.lap[3,effsize[[3]],2],
            theta.lap[4,effsize[[4]],2])
)

data.lap1$group <- factor(data.lap1$group,
                          levels = c("24800","37750","53700","92800")
)
beta1lapbox <- ggplot(data.lap1, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 1, color ="red") + 
  scale_y_continuous(limits = c(0.9, 1.1)) +
  labs(title=expression(beta[2]~"glmmTMB"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.regva1 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.regva[1,effsize[[1]],3],
            theta.regva[2,effsize[[2]],3],
            theta.regva[3,effsize[[3]],3],
            theta.regva[4,effsize[[4]],3])
)

data.regva1$group <- factor(data.regva1$group,
                            levels = c("24800","37750","53700","92800")
)
beta1regvabox <- ggplot(data.regva1, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 1, color ="red") + 
  scale_y_continuous(limits = c(0.9,1.1)) +
  labs(title=expression(beta[2]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

library(patchwork)
beta0gvabox  + beta0lapbox +beta0regvabox +
  beta1gvabox + beta1lapbox +beta1regvabox + plot_layout(ncol = 3)

data.gva2 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.gva[1,effsize[[1]],3],
            theta.gva[2,effsize[[2]],3],
            theta.gva[3,effsize[[3]],3],
            theta.gva[4,effsize[[4]],3])
)

data.gva2$group <- factor(data.gva2$group,
                          levels = c("24800","37750","53700","92800")
)
beta2gvabox <- ggplot(data.gva2, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -0.5, color ="red") + 
  scale_y_continuous(limits = c(-0.6,-0.4)) +
  labs(title=expression(beta[3]~"GVA"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.lap2 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.lap[1,effsize[[1]],3],
            theta.lap[2,effsize[[2]],3],
            theta.lap[3,effsize[[3]],3],
            theta.lap[4,effsize[[4]],3])
)

data.lap2$group <- factor(data.lap2$group,
                          levels = c("24800","37750","53700","92800")
)
beta2lapbox <- ggplot(data.lap2, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -0.5, color ="red") + 
  scale_y_continuous(limits = c(-0.6,-0.4)) +
  labs(title=expression(beta[3]~"glmmTMB"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.regva2 <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.regva[1,effsize[[1]],4],
            theta.regva[2,effsize[[2]],4],
            theta.regva[3,effsize[[3]],4],
            theta.regva[4,effsize[[4]],4])
)

data.regva2$group <- factor(data.regva2$group,
                            levels = c("24800","37750","53700","92800")
)
beta2regvabox <- ggplot(data.regva2, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -0.5, color ="red") + 
  scale_y_continuous(limits = c(-0.6,-0.4)) +
  labs(title=expression(beta[3]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

###boxplot for beta###
library(patchwork)
beta0gvabox  + beta0regvabox + beta0lapbox +  
  beta1gvabox + beta1regvabox + beta1lapbox +
  beta2gvabox + beta2regvabox + beta2lapbox  + plot_layout(ncol = 3)

data.gvau <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.gva[1,effsize[[1]],p+1],
            theta.gva[2,effsize[[2]],p+1],
            theta.gva[3,effsize[[3]],p+1],
            theta.gva[4,effsize[[4]],p+1])
)
data.gvau$group <- factor(data.gvau$group,
                          levels = c("24800","37750","53700","92800")
)
sigma1gvabox <- ggplot(data.gvau, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") +
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.25, 0.75)) +
  labs(title=expression(sigma[u]~"GVA"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.lapu <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.lap[1,effsize[[1]],p+1],
            theta.lap[2,effsize[[2]],p+1],
            theta.lap[3,effsize[[3]],p+1],
            theta.lap[4,effsize[[4]],p+1])
)
data.lapu$group <- factor(data.lapu$group,
                          levels = c("24800","37750","53700","92800")
)
sigma1lapbox <- ggplot(data.lapu, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.25, 0.75)) +
  labs(title=expression(sigma[u]~"glmmTMB"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.regvau <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.regva[1,effsize[[1]],p+2],
            theta.regva[2,effsize[[2]],p+2],
            theta.regva[3,effsize[[3]],p+2],
            theta.regva[4,effsize[[4]],p+2])
)

data.regvau$group <- factor(data.regvau$group,
                            levels = c("24800","37750","53700","92800")
)
sigma1regvabox <- ggplot(data.regvau, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.25, 0.75)) +
  labs(title=expression(sigma[u]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.gvav <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.gva[1,effsize[[1]],p+2],
            theta.gva[2,effsize[[2]],p+2],
            theta.gva[3,effsize[[3]],p+2],
            theta.gva[4,effsize[[4]],p+2])
)

data.gvav$group <- factor(data.gvav$group,
                          levels = c("24800","37750","53700","92800")
)
sigma2gvabox <- ggplot(data.gvav, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.3, 0.7)) +
  labs(title=expression(sigma[v]~"GVA"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

data.lapv <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.lap[1,effsize[[1]],p+2],
            theta.lap[2,effsize[[2]],p+2],
            theta.lap[3,effsize[[3]],p+2],
            theta.lap[4,effsize[[4]],p+2])
)

data.lapv$group <- factor(data.lapv$group,
                          levels = c("24800","37750","53700","92800")
)
sigma2lapbox <- ggplot(data.lapv, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.3, 0.7)) +
  labs(title=expression(sigma[v]~"glmmTMB"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))


data.regvav <- data.frame(
  group = c(rep("24800", length(effsize[[1]])),
            rep("37750", length(effsize[[2]])), 
            rep("53700", length(effsize[[3]])),
            rep("92800", length(effsize[[4]]))),
  value = c(theta.regva[1,effsize[[1]],p+3],
            theta.regva[2,effsize[[2]],p+3],
            theta.regva[3,effsize[[3]],p+3],
            theta.regva[4,effsize[[4]],p+3])
)

data.regvav$group <- factor(data.regvav$group,
                            levels = c("24800","37750","53700","92800")
)
sigma2regvabox <- ggplot(data.regvav, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") +
  scale_y_continuous(limits = c(0.3, 0.7)) +
  labs(title=expression(sigma[v]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1))

library(patchwork)
sigma1gvabox  + sigma1regvabox + sigma1lapbox +
  sigma2gvabox + sigma2regvabox + sigma2lapbox + plot_layout(ncol = 3)



