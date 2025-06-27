########## Simulation for Gamma case ##########
# Last changed: 20 January. 2025

###############################################################################
rm(list = ls())
##GAMMA case##
#setwd("xx\\code\\Gamma") ###set the path###
setwd("C:\\Users\\DELL\\Desktop\\Resubmit-SS-R3\\Gamma")
library(Matrix)
library(MASS)
library(statmod)
library(mvQuad)
library(mvtnorm)

###############################################################################
source("CalculateB.Rs")
source("CLGVA.R")
source("GVAuv.R")
do.GVAUiVj <- TRUE
do.GVACRCL <- TRUE

###############################################################################
###Simulation settings###
alpha <- 0.8
Sigma <- c(0.5,0.5)
delta <- 1.1
M <- c(18,20,30,50,80,100,150,200,300)
N <- round(M^(1/delta),0)

beta.true <-  c(-2,1,-0.5,1,-1,0.5,0.8,-0.4,0.1)
p <- length(beta.true)

TRIALS  <- 500  # Total number of trials
FAMILY  <- c("GAMMA")
###############################################################################
###Record results
theta.gva          <- matrix(0,length(M)*TRIALS*(p+2),1); dim(theta.gva)          <- c(length(M),TRIALS,p+2)
theta.regva        <- matrix(0,length(M)*TRIALS*(p+3),1); dim(theta.regva)        <- c(length(M),TRIALS,p+3)

se.theta.gva       <- matrix(0,length(M)*TRIALS*(p+2),1); dim(se.theta.gva)       <- c(length(M),TRIALS,p+2)
se.theta.regva     <- matrix(0,length(M)*TRIALS*(p+3),1); dim(se.theta.regva)     <- c(length(M),TRIALS,p+3)
asy.se.theta.regva <- matrix(0,length(M)*TRIALS*(p+2),1); dim(asy.se.theta.regva) <- c(length(M),TRIALS,p+2)

times.gva    <- matrix(0,length(M)*TRIALS,1); dim(times.gva) <- c(length(M),TRIALS)
times.regva  <- matrix(0,length(M)*TRIALS,1); dim(times.regva) <- c(length(M),TRIALS)

converged.gva    <- matrix(0,length(M)*TRIALS); dim(converged.gva) <- c(length(M),TRIALS)
converged.regva  <- matrix(0,length(M)*TRIALS); dim(converged.regva) <- c(length(M),TRIALS)

################################################
ITER0 <- 1
trial0 <- 1
for (ITER in ITER0:length(M)) {
  for (trial in trial0:TRIALS) {
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
    
    vmu.true <- exp(veta)
    scale0 <- vmu.true/alpha
    rp <- rgamma(n*m,alpha,scale = scale0)
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
    
    family = "GAMMA"
    if (do.GVACRCL){
      cat("trial=",trial,"    m=",m," n=",n,"    sigmaw=",Sigma[1]," sigmav=",Sigma[2],"\n")
      cat("Fitting model using ReGVA.CLCR method \n")
      
      bval1  <- proc.time()
      res <- try(ReGVA.CLCR.fit <- ReGVA.CLCR.FIT(alpha,vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id),silent=TRUE)
      eval1  <- proc.time()
      
      if ((class(res)[1])!="try-error"){
        theta.regva[ITER,trial,1:(p+1)] <- res$vbeta
        theta.regva[ITER,trial,p+2:3] <- sqrt(res$sigma2)
        times.regva[ITER,trial]  <- eval1[3]-bval1[3]	
        ##calculate the se by aym...
        sumx <- matrix(0,m,n) 
        XX <- matrix(0,m*n,p-1)
        for(i in 1:(p-1)){
          XX[,i] <- c(mX[,,i])
        }
        
        SigmaXinv <- solve(cov(XX))
        sesigma1 <- sqrt(res$sigma2[1]/(2*m))
        sesigma2 <- sqrt(res$sigma2[2]/(2*n))
        
        a11 <- res$sigma2[2]/4
        a12 <- res$sigma2[2]/4
        a13 <- 0
        a22 <- (exp(res$sigma2[2])-1)/4
        a23 <- -res$sigma2[2]^2/8
        a33 <- res$sigma2[2]^2/8
        Gamma1 <- matrix(c(a11,a12,a13,a12,a22,a23,a13,a23,a33),3,3) 
        b11 <- res$sigma2[1]/4
        b21 <- res$sigma2[1]/4
        b31 <- 0
        b22 <- (exp(res$sigma2[1])-1)/4
        b23 <- -res$sigma2[1]^2/8
        b33 <- res$sigma2[1]^2/8
        Gamma2 <- matrix(c(b11,b21,b31,b21,b22,b23,b31,b23,b33),3,3) 
        
        I <- c(1,1,1)
        sebeta0 <- sqrt(t(I)%*%Gamma1%*%I/(n)+t(I)%*%Gamma2%*%I/(m))
        
        c11 <- ((exp(res$sigma2[2])-1) + exp(res$sigma2[2])/alpha)/4
        c12 <- 1/(4*alpha)
        c21 <- 1/(4*alpha)
        c22 <- ((exp(res$sigma2[1])-1) + exp(res$sigma2[1])/alpha)/4
        
        Gamma3 <- matrix(c(c11,c21,c12,c22),2,2)
        I2 <- c(1,1)
        sebeta1 <- sqrt(diag(SigmaXinv))*c(sqrt(t(I2)%*%Gamma3%*%I2/(n*m)))
        
        se.theta.regva[ITER,trial,]     <- res$vbeta.serr
        asy.se.theta.regva[ITER,trial,] <- c(sebeta0,sebeta1,sesigma1,sesigma2) 
        converged.regva[ITER,trial] <- TRUE    
      } else {
        converged.regva[ITER,trial] <- FALSE
      }
    }
    
    ###############################################################################
    if (do.GVAUiVj){
      cat("trial=",trial,"    m=",m,"    sigmaw=",Sigma[1]," sigmav=",Sigma[2],"\n")
      cat("Fitting model using GVA.UiVj \n")
      
      vy <- t(matrix(rp,n,m))
      mX <- X
      #############################################################################
      vbeta  <- rep(0,p)
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
      res1 <- try(GVA.UiVj.fit <- GVA.UiVj.FIT(alpha,vbeta,vSigma,lumu,luLambda,lvmu,lvLambda,family,vy,mX,id),silent=TRUE)
      eval2  <- proc.time()
      
      if ((class(res)[1])!="try-error") {			
        theta.gva[ITER,trial,1:p] <- res1$vbeta
        theta.gva[ITER,trial,p+(1:2)] <- sqrt(res1$sigma2)
        
        se.theta.gva[ITER,trial,1:(p+2)] <- res1$vbeta.serr
        times.gva[ITER,trial] <- eval2[3]-bval2[3]	
        converged.gva[ITER,trial] <- TRUE    
      } else {
        converged.gva[ITER,trial] <- FALSE
      }	
    }
    ###############################################################################
    save.image("Gam-S-M500-delta1-1-1.Rdata")	
  }
}

load("C:\\Users\\DELL\\Desktop\\Resubmit-SS-R3\\Gamma\\Gam-S-M500-delta1-1-1.Rdata")
###############################################################################
indOutliers <- function(x,nsd=5) {
  z.scores <- (x - median(x))/mad(x)
  return(which(abs(z.scores)>nsd))
}

gva.beta0 <- list()
gva.beta1 <- list()
gva.sigma1 <- list()
gva.sigma2 <- list()

regva.beta0 <- list()
regva.beta1 <- list()
regva.sigma1 <- list()
regva.sigma2 <- list()

gva.time <- c()
regva.time <- c()

effsize <- list()


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
  
  inds.outliers <- c()
  for(i in 1:dim(theta.gva)[3]){
    inds.outliers <- c(inds.outliers,
                       indOutliers(theta.regva[ITER,,i]),
                       indOutliers(theta.gva[ITER,,i]),
                       which(asy.se.theta.regva[ITER,,i]== "NaN"),
                       which(se.theta.gva[ITER,,i]=="NaN")
    )
  }
  inds.outliers <- c(inds.outliers,
                     indOutliers(theta.regva[ITER,,p+3]),
                     which(se.theta.regva[ITER,,p+3]=="NaN"))
  
  inds.outliers <- c(uncong.gva,uncong.regva,inds.outliers)
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
  mean.beta01.regva <- round(locfun(theta.regva[ITER,inds,1]),dp)
  mean.beta02.regva <- round(locfun(theta.regva[ITER,inds,2]),dp)
  mean.beta0.regva  <- round(locfun((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,(p+2)]^2/2-theta.regva[ITER,inds,(p+3)]^2/2)/2),dp)
  
  meseasy.beta0.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,1]),dp)
  
  sd.beta0.gva     <- round(sprfun(theta.gva[ITER,inds,1]),dp)  
  sd.beta0.regva   <- round(sprfun((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,(p+2)]^2/2-theta.regva[ITER,inds,(p+3)]^2/2)/2),dp)
  
  mstderr.beta0.gva    <- round(locfun(se.theta.gva[ITER,inds,1]),dp)
  
  rmse.beta0.gva <-  round(sqrt(mean((theta.gva[ITER,inds,1] - beta.true[1])^{2})),dp)
  rmse.beta0.regva <- round(sqrt(locfun(((theta.regva[ITER,inds,1]+theta.regva[ITER,inds,2]-theta.regva[ITER,inds,(p+2)]^2/2-theta.regva[ITER,inds,(p+3)]^2/2)/2 - beta.true[1])^2)),dp)
  
  cat("beta0.gva \n")
  
  cat(mean.beta0.gva, "(" ,sd.beta0.gva , ") & " , mstderr.beta0.gva, "&", rmse.beta0.gva, "\n") 
  cat("beta0.regva \n")
  cat(mean.beta0.regva, "(",sd.beta0.regva ,") & ", meseasy.beta0.regvaasy ,"&", rmse.beta0.regva, " \n")  
  
  cat("beta1 \n")
  
  mean.beta1.gva   <- round(locfun(theta.gva[ITER,inds,2]),dp)
  mean.beta1.regva <- round(locfun(theta.regva[ITER,inds,3]),dp)
  
  sd.beta1.gva   <- round(sprfun(theta.gva[ITER,inds,2]),dp)  
  sd.beta1.regva <- round(sprfun(theta.regva[ITER,inds,3]),dp)  
  
  mstderr.beta1.gva   <- round(locfun(se.theta.gva[ITER,inds,2]),dp)
  meseasy.beta1.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,2]),dp)
  
  rmse.beta1.gva <-  round(sqrt(mean((theta.gva[ITER,inds,2] - beta.true[2])^{2})),dp)
  rmse.beta1.regva <- round(sqrt(locfun((theta.regva[ITER,inds,3] - beta.true[2])^{2})),dp)
  
  
  cat("beta1.gva \n")
  cat(mean.beta1.gva,"(",sd.beta1.gva ,") & ", mstderr.beta1.gva, "& ", rmse.beta1.gva, " \n")
  cat("beta1.regva \n")
  cat(mean.beta1.regva, "(",sd.beta1.regva ,") & ", meseasy.beta1.regvaasy, "& ", rmse.beta1.regva, " \n")
  
  cat("sigma \n")
  mean.sigmau.gva   <- round(locfun(theta.gva[ITER,inds,p+1]),dp)
  mean.sigmav.gva   <- round(locfun(theta.gva[ITER,inds,p+2]),dp)
  
  mean.sigmau.regva <- round(locfun(theta.regva[ITER,inds,p+2]),dp)
  mean.sigmav.regva <- round(locfun(theta.regva[ITER,inds,p+3]),dp)
  
  sd.sigmau.gva   <- round(sprfun(theta.gva[ITER,inds,p+1]),dp)  
  sd.sigmav.gva   <- round(sprfun(theta.gva[ITER,inds,p+2]),dp)
  
  sd.sigmau.regva <- round(sprfun(theta.regva[ITER,inds,p+2]),dp)
  sd.sigmav.regva <- round(sprfun(theta.regva[ITER,inds,p+3]),dp)
  
  
  mstderr.sigma1.gva <-  round(locfun(se.theta.gva[ITER,inds,p+1]*theta.gva[ITER,inds,p+1]/2),dp)
  mstderr.sigma2.gva <-  round(locfun(se.theta.gva[ITER,inds,p+2]*theta.gva[ITER,inds,p+2]/2),dp)
  
  meseasy.sigma1.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,p+1]),dp)
  meseasy.sigma2.regvaasy <- round(locfun(asy.se.theta.regva[ITER,inds,p+2]),dp)
  
  rmse.sigma1.gva <-  round(sqrt(mean((theta.gva[ITER,inds,p+1] - Sigma[1])^{2})),dp)
  rmse.sigma2.gva <-  round(sqrt(mean((theta.gva[ITER,inds,p+2] - Sigma[2])^{2})),dp)
  
  rmse.sigma1.regva <- round(sqrt(mean((theta.regva[ITER,inds,p+2] - Sigma[1])^{2})),dp)
  rmse.sigma2.regva <- round(sqrt(mean((theta.regva[ITER,inds,p+3] - Sigma[2])^{2})),dp)
  
  cat(mean.sigmau.gva , "(",sd.sigmau.gva ,") &  ", mstderr.sigma1.gva, " & ", rmse.sigma1.gva, " \n") 
  cat(mean.sigmav.gva , "(",sd.sigmav.gva ,") &  ", mstderr.sigma2.gva, " & ", rmse.sigma2.gva, " \n") 
  
  cat(mean.sigmau.regva ,"(",sd.sigmau.regva ,") &  ", meseasy.sigma1.regvaasy, " & ", rmse.sigma1.regva," \n")
  cat(mean.sigmav.regva ,"(",sd.sigmav.regva ,") &  ", meseasy.sigma2.regvaasy, " & ", rmse.sigma2.regva," \n")
  
  cat("times - ", round(mean(times.gva[ITER,inds]),dp),
      round(mean(times.regva[ITER,inds]),dp),"\n") 
  
  cat("\n\n") 
  
  gva.time <- c(gva.time, round(mean(times.gva[ITER,inds]),dp))
  regva.time <- c(regva.time, round(mean(times.regva[ITER,inds]),dp))
}



###plot the figure of CPU Time vs sample size###

###plot Figure 1 of CPU Time vs sample size###
gvatime <- data.frame(x=M*N,y=gva.time)
regvatime <- data.frame(x=M*N,y=regva.time)
library(ggplot2)
log10_labels <- function(x) {
  parse(text = paste0("10^", log10(x)))  # 转换为 10^x 格式
}
plot1 <- 
  ggplot()+                   # 添加连线
  geom_point(data = gvatime, aes(x=x, y=y, color = "blue"), size =3) + 
  geom_point(data = regvatime, aes(x=x, y=y, color = "green"), size =3) +
  geom_line(data = gvatime, aes(x=x, y=y, color = "blue")) +
  geom_line(data = regvatime, aes(x=x, y=y, color = "green")) +
  scale_color_manual(values =c("blue","green"),labels = c("GVA","GVACL")) +
  scale_y_continuous(
    trans = "log10",                # 纵轴对数变换
    limits = c(4e-2, 5e2),          # 设置纵轴范围
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
    title = "CPU sec. Gamma",     # 图表标题
    x = "Sample size",              # 横轴标签
    y = "Time/s",                   # 纵轴标签
    color = "Method"                # 图例标题
  )+
  labs(title = "CPU sec. Gamma", x= "Sample size", y="Time/s") +
  theme_minimal() + 
  theme(legend.position = "top")


gvatime <- data.frame(x=log(M*N,10),y=log(gva.time,10))
regvatime <- data.frame(x=log(M*N,10),y=log(regva.time,10))
library(ggplot2)
plot1 <- ggplot() + 
  geom_point(data = gvatime, aes(x=x, y=y, color = "blue"), size =3) + 
  geom_point(data = regvatime, aes(x=x, y=y, color = "green"), size =3) +
  geom_line(data = gvatime, aes(x=x, y=y, color = "blue")) +
  geom_line(data = regvatime, aes(x=x, y=y, color = "green")) +
  scale_color_manual(values =c("blue","green"),labels = c("GVA","GVACL")) +
  scale_x_continuous(breaks = log10(c(1e2,1e3,1e4,1e5)),
                     labels = scales::label_math(10^.x)) +
  scale_y_continuous(breaks = log10(c(1e-2,1e-1,1e0,1e1,1e2)),
                     labels = scales::label_math(10^.x)) +
  labs(title = "CPU sec.  Gamma", x= "Sample size", y="Time/s") +
  theme_minimal() + 
  theme(legend.position = "top")

###plot boxplots of parameter estimates###
data.gva0 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.gva[1,effsize[[1]],1], 
            theta.gva[2,effsize[[2]],1],
            theta.gva[3,effsize[[3]],1],
            theta.gva[4,effsize[[4]],1],
            theta.gva[5,effsize[[5]],1],
            theta.gva[6,effsize[[6]],1],
            theta.gva[7,effsize[[7]],1],
            theta.gva[8,effsize[[8]],1],
            theta.gva[9,effsize[[9]],1])
)

data.gva0$group <- factor(data.gva0$group,
                          levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta0gvabox <- ggplot(data.gva0, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") +
  geom_hline( yintercept = -2, color ="red") + 
  scale_y_continuous(limits = c(-2.5, -1.5)) +
  labs(title=expression(beta[1]~"GVA"), x="Sample size ", y="" ) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.regva0 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c((theta.regva[1,effsize[[1]],1] + theta.regva[1,effsize[[1]],2] - theta.regva[1,effsize[[1]],p+2]^2/2-theta.regva[1,effsize[[1]],p+3]^2/2)/2,
            (theta.regva[2,effsize[[2]],1] + theta.regva[2,effsize[[2]],2] - theta.regva[2,effsize[[2]],p+2]^2/2-theta.regva[2,effsize[[2]],p+3]^2/2)/2,
            (theta.regva[3,effsize[[3]],1] + theta.regva[3,effsize[[3]],2] - theta.regva[3,effsize[[3]],p+2]^2/2-theta.regva[3,effsize[[3]],p+3]^2/2)/2,
            (theta.regva[4,effsize[[4]],1] + theta.regva[4,effsize[[4]],2] - theta.regva[4,effsize[[4]],p+2]^2/2-theta.regva[4,effsize[[4]],p+3]^2/2)/2,
            (theta.regva[5,effsize[[5]],1] + theta.regva[5,effsize[[5]],2] - theta.regva[5,effsize[[5]],p+2]^2/2-theta.regva[5,effsize[[5]],p+3]^2/2)/2,
            (theta.regva[6,effsize[[6]],1] + theta.regva[6,effsize[[6]],2] - theta.regva[6,effsize[[6]],p+2]^2/2-theta.regva[6,effsize[[6]],p+3]^2/2)/2,
            (theta.regva[7,effsize[[7]],1] + theta.regva[7,effsize[[7]],2] - theta.regva[7,effsize[[7]],p+2]^2/2-theta.regva[7,effsize[[7]],p+3]^2/2)/2,
            (theta.regva[8,effsize[[8]],1] + theta.regva[8,effsize[[8]],2] - theta.regva[8,effsize[[8]],p+2]^2/2-theta.regva[8,effsize[[8]],p+3]^2/2)/2,
            (theta.regva[9,effsize[[9]],1] + theta.regva[9,effsize[[9]],2] - theta.regva[9,effsize[[9]],p+2]^2/2-theta.regva[9,effsize[[9]],p+3]^2/2)/2)
)

data.regva0$group <- factor(data.regva0$group,
                            levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta0regvabox <- ggplot(data.regva0, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -2, color ="red") + 
  scale_y_continuous(limits = c(-2.5, -1.5)) +
  labs(title=expression(beta[1]~"GVACL"), x="Sample size", y="") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )
beta0gvabox + beta0regvabox

data.gva1 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.gva[1,effsize[[1]],2], 
            theta.gva[2,effsize[[2]],2],
            theta.gva[3,effsize[[3]],2],
            theta.gva[4,effsize[[4]],2],
            theta.gva[5,effsize[[5]],2],
            theta.gva[6,effsize[[6]],2],
            theta.gva[7,effsize[[7]],2],
            theta.gva[8,effsize[[8]],2],
            theta.gva[9,effsize[[9]],2])
)

data.gva1$group <- factor(data.gva1$group,
                          levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta1gvabox <- ggplot(data.gva1, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 1, color ="red") +
  scale_y_continuous(limits = c(0.8, 1.2)) +
  labs(title=expression(beta[2]~"GVA"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )


data.regva1 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.regva[1,effsize[[1]],3],
            theta.regva[2,effsize[[2]],3],
            theta.regva[3,effsize[[3]],3],
            theta.regva[4,effsize[[4]],3],
            theta.regva[5,effsize[[5]],3],
            theta.regva[6,effsize[[6]],3],
            theta.regva[7,effsize[[7]],3],
            theta.regva[8,effsize[[8]],3],
            theta.regva[9,effsize[[9]],3])
)

data.regva1$group <- factor(data.regva1$group,
                            levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta1regvabox <- ggplot(data.regva1, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 1, color ="red") + 
  scale_y_continuous(limits = c(0.8, 1.2)) +
  labs(title=expression(beta[2]~"GVACL"), x="Sample size ", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.gva2 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.gva[1,effsize[[1]],3], 
            theta.gva[2,effsize[[2]],3],
            theta.gva[3,effsize[[3]],3],
            theta.gva[4,effsize[[4]],3],
            theta.gva[5,effsize[[5]],3],
            theta.gva[6,effsize[[6]],3],
            theta.gva[7,effsize[[7]],3],
            theta.gva[8,effsize[[8]],3],
            theta.gva[9,effsize[[9]],3])
)

data.gva2$group <- factor(data.gva2$group,
                          levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta2gvabox <- ggplot(data.gva2, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") +
  geom_hline( yintercept = -0.5, color ="red") + 
  scale_y_continuous(limits = c(-0.7, -0.3)) +
  labs(title=expression(beta[3]~"GVA"), x="Sample size", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.regva2 <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.regva[1,effsize[[1]],4],
            theta.regva[2,effsize[[2]],4],
            theta.regva[3,effsize[[3]],4],
            theta.regva[4,effsize[[4]],4],
            theta.regva[5,effsize[[5]],4],
            theta.regva[6,effsize[[6]],4],
            theta.regva[7,effsize[[7]],4],
            theta.regva[8,effsize[[8]],4],
            theta.regva[9,effsize[[9]],4])
)

data.regva2$group <- factor(data.regva2$group,
                            levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
beta2regvabox <- ggplot(data.regva2, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = -0.5, color ="red") + 
  scale_y_continuous(limits = c(-0.7, -0.3)) +
  labs(title=expression(beta[3]~"GVACL"), x="Sample size", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

###boxplot of beta##
library(patchwork)
beta0gvabox  + beta0regvabox +
  beta1gvabox + beta1regvabox +
  beta2gvabox + beta2regvabox + plot_layout(ncol = 2)


data.gvau <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.gva[1,effsize[[1]],p+1],
            theta.gva[2,effsize[[2]],p+1],
            theta.gva[3,effsize[[3]],p+1],
            theta.gva[4,effsize[[4]],p+1],
            theta.gva[5,effsize[[5]],p+1],
            theta.gva[6,effsize[[6]],p+1],
            theta.gva[7,effsize[[7]],p+1],
            theta.gva[8,effsize[[8]],p+1],
            theta.gva[9,effsize[[9]],p+1])
)
data.gvau$group <- factor(data.gvau$group,
                          levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
sigma1gvabox <- ggplot(data.gvau, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") +
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.1, 1)) +
  labs(title=expression(sigma[u]~"GVA"), x="Sample size", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )


data.regvau <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.regva[1,effsize[[1]],p+2],
            theta.regva[2,effsize[[2]],p+2],
            theta.regva[3,effsize[[3]],p+2],
            theta.regva[4,effsize[[4]],p+2],
            theta.regva[5,effsize[[5]],p+2],
            theta.regva[6,effsize[[6]],p+2],
            theta.regva[7,effsize[[7]],p+2],
            theta.regva[8,effsize[[8]],p+2],
            theta.regva[9,effsize[[9]],p+2])
)

data.regvau$group <- factor(data.regvau$group,
                            levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
sigma1regvabox <- ggplot(data.regvau, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.1, 1)) +
  labs(title=expression(sigma[u]~"GVACL"), x="Sample size", y= "") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.gvav <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600",  length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.gva[1,effsize[[1]],p+2],
            theta.gva[2,effsize[[2]],p+2],
            theta.gva[3,effsize[[3]],p+2],
            theta.gva[4,effsize[[4]],p+2],
            theta.gva[5,effsize[[5]],p+2],
            theta.gva[6,effsize[[6]],p+2],
            theta.gva[7,effsize[[7]],p+2],
            theta.gva[8,effsize[[8]],p+2],
            theta.gva[9,effsize[[9]],p+2])
)

data.gvav$group <- factor(data.gvav$group,
                          levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
sigma2gvabox <- ggplot(data.gvav, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.1, 1)) +
  labs(title=expression(sigma[v]~"GVA"), x="Sample size", y="" ) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

data.regvav <- data.frame(
  group = c(rep("252",   length(effsize[[1]])),
            rep("300",   length(effsize[[2]])),
            rep("660",   length(effsize[[3]])),
            rep("1750",  length(effsize[[4]])),
            rep("4320",  length(effsize[[5]])),
            rep("6600", length(effsize[[6]])),
            rep("14250", length(effsize[[7]])), 
            rep("24800", length(effsize[[8]])),
            rep("53700", length(effsize[[9]]))),
  value = c(theta.regva[1,effsize[[1]],p+3],
            theta.regva[2,effsize[[2]],p+3],
            theta.regva[3,effsize[[3]],p+3],
            theta.regva[4,effsize[[4]],p+3],
            theta.regva[5,effsize[[5]],p+3],
            theta.regva[6,effsize[[6]],p+3],
            theta.regva[7,effsize[[7]],p+3],
            theta.regva[8,effsize[[8]],p+3],
            theta.regva[9,effsize[[9]],p+3])
)

data.regvav$group <- factor(data.regvav$group,
                            levels = c("252","300","660","1750","4320","6600","14250","24800","53700")
)
sigma2regvabox <- ggplot(data.regvav, aes(x=group, y=value, fill = group)) + 
  geom_boxplot(fill="gray", outlier.color = "gray", outlier.fill="gray") + 
  geom_hline( yintercept = 0.5, color ="red") + 
  scale_y_continuous(limits = c(0.1, 1)) +
  labs(title=expression(sigma[v]~"GVACL"), x="Sample size", y="" ) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

###boxplot of sigma###
library(patchwork)
sigma1gvabox  + sigma1regvabox +
  sigma2gvabox + sigma2regvabox + plot_layout(ncol = 2)



d <- dim(asy.se.theta.regva)[3]
CP <- matrix(0,length(M),d) 
para <- c(beta.true,Sigma)
for(j in 1:p){
  for(i in 1:length(M)){
    #asy.se.theta.regva[i,1:trial,1]
    #dim(theta.gva)
    UCI <- theta.gva[i,1:trial,j] + 1.96*asy.se.theta.regva[i,1:trial,j]
    LCI <- theta.gva[i,1:trial,j] - 1.96*asy.se.theta.regva[i,1:trial,j]
    UI <- which(UCI>=para[j])
    LI <- which(LCI<=para[j])
    CP[i,j] <- 1 - length(setdiff(UI,LI))/trial
  }
}
# 加载必要的库
library(ggplot2)
M*N
# 示例数据（你需要用实际数据替换这里的假数据）
set.seed(123)
m_values <- M*N
beta0_coverage <- CP[,1]*100  # 模拟β0的覆盖率
beta1_coverage <- CP[,2]*100  # 模拟β1的覆盖率
beta2_coverage <- CP[,3]*100
sigma1_sq_coverage <- CP[,p-1]*100  # 模拟σ^2的覆盖率
sigma2_sq_coverage <- CP[,p]*100  # 模拟σ^2的覆盖率

data <- data.frame(
  m = rep(m_values, 5),
  coverage = c(beta0_coverage, beta1_coverage, beta2_coverage, sigma1_sq_coverage,sigma2_sq_coverage),
  parameter = factor(rep(c("beta0", "beta1", "beta2","sigmau","sigmav"), each = length(m_values)))
)

#data <- data.frame(
#  m = rep(m_values, 5),
#  coverage = c(beta0_coverage, beta1_coverage, beta2_coverage, sigma1_sq_coverage,sigma2_sq_coverage),
#  parameter = factor(rep(c(beta[0], beta[1], beta[2],sigma[u]^2,sigma[v]^2), each = length(m_values)))
#)

# 创建图形
GamCPplot <- ggplot(data, aes(x = m, y = coverage, linetype = parameter)) +
  geom_line() +
  geom_hline(yintercept = 95, linetype = "solid", color = "grey", size = 1) +  # 加粗参考线
  labs(
    x = "sample size",
    y = "coverage percentage",
    title = "Gamma case"~~~~X[ij] %~% N(0, 1),
    subtitle = expression(beta[0]^0 == -2 ~ "," ~ beta[1]^0 == 1 ~ ","~beta[2]^0==-0.5~"," 
                          ~ sigma[u]^0 == 0.5~","~ sigma[v]^0 == 0.5),
    linetype = "Parameter"  # 图例标题
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted","dotdash","longdash"),
    labels = c(expression(beta[0]^0), expression(beta[1]^0), 
               expression(beta[2]^0), expression(sigma[u]^0),
               expression(sigma[v]^0))  # 使用希腊字母作为标签
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  ylim(85, 100) +  # 设置y轴最大值为100
  scale_x_continuous(
    trans = "log10",                # 横轴对数变换
    limits = c(1e2, 1e5),           # 设置横轴范围
    #breaks = c(1e2, 1e3), # 设置刻度
    labels = log10_labels           # 使用 10^x 格式标签
  )
GamCPplot

library(patchwork)
PoiCPplot  + GamCPplot +  plot_layout(ncol = 2)

