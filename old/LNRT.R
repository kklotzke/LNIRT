
LNRT <- function(RT,XG,Discrimination,WL){

  N <- nrow(RT) 
  K <- ncol(RT)  

  #time discrimination
  if(missing(Discrimination)){
	Discrimination <- TRUE
  }else{
	if(Discrimination == TRUE){
		Discrimination <- TRUE
	}else{
		Discrimination <- FALSE	
	}
  } 

  if (missing(WL)) {
	WL <- 0 
  }else{
 	WL <- 1	#time discrimination = 1/sqrt(error variance)
  }

  ## population theta (ability - speed)
  theta <- matrix(rnorm(N))
  muP <- matrix(0,N,1)
  SigmaP <-  diag(1)
  muP0 <- matrix(0,1,1)
  SigmaP0 <-  diag(1)/100
  ingroup <- rep(1,N) 
  Mingroup <-  matrix(0,ncol=2,nrow=N)
  
  if(XG > 1000){    
    flagged <- matrix(0,ncol=1,nrow=N)     
    ss <- 1
  }  
  
  ##population item (ability - speed)
  ab <- matrix(rnorm(K*2),ncol=2)
  ab[,1] <- 1
  muI <- t(matrix(rep(c(1,0),K),ncol=K))
  muI0 <- muI[1,]
  SigmaI <- diag(2)
  SigmaI0 <- diag(2)*10
  for(ii in 1:2) {
    SigmaI0[ii,] <- SigmaI0[ii,]*c(.5,3)
  }
  
  ##storage
  MT <- MT2 <- array(0,dim=c(N))
  MAB <- array(0,dim=c(XG,K,2))
  MmuP <- array(0,dim=c(XG,1))
  MmuI <- array(0,dim=c(XG,2))
  MSP <- array(0,dim=c(XG,1,1))
  MSI <- array(0,dim=c(XG,2,2))
  sigma2 <- rep(1,K)
  Msigma2 <- matrix(0,ncol=K,nrow=XG)
  lZP <- 0
  lZPT <- 0
  lZI <- 0	
  EAPresid <- matrix(0,ncol=K,nrow=N)
  EAPKS <- matrix(0,ncol=1,nrow=K)
  iis <- 1
  EAPphi <- matrix(0,ncol=1,nrow=K)
  EAPlambda <- matrix(0,ncol=1,nrow=K)
  EAPtheta <- matrix(0,ncol=1,nrow=N)
  EAPsigma2	<- matrix(0,ncol=1,nrow=K)

  DT 	<- matrix(1,ncol=1,nrow=N*K) 
  DT[which(is.na(RT))] <- 0 
  DT <- matrix(DT,nrow=N,ncol=K)

  EAPCP <- matrix(0,ncol=1,nrow=N)	
  
  ## Start MCMC algorithm
  
  for (ii in 1:XG){
	
	if(sum(DT==0) > 0){
		if(WL == 1){
			RT <- SimulateRT(RT=RT,zeta=theta,lambda=ab[,2],phi=rep(1,K),sigma2=sigma2,DT=DT)
		}else{
			RT <- SimulateRT(RT=RT,zeta=theta,lambda=ab[,2],phi=ab[,1],sigma2=sigma2,DT=DT)
		}
	}
    
    theta <- DrawZeta_LNRT(RT,ab[,1],ab[,2],sigma2,muP,SigmaP[1,1]) 
    theta[1:N] <- theta[1:N] - mean(theta) 
    MT[1:N] <- MT[1:N] + theta[1:N]
    MT2[1:N] <- MT2[1:N] + theta[1:N]^2 
   
	if((WL == 1)){
		#no time discrimination, 1/(sqrt(error variance)) = discrimination on MVN prior 
		ab[,1] <- 1	
		ab1 <- cbind(1/sqrt(sigma2),ab[,2])
    		dum <- Conditional(kk=2,Mu=muI,Sigma=SigmaI,Z=ab1)
		ab[,2] <- DrawLambda(RT=RT,zeta=theta,sigma2=sigma2,mu=dum$CMU,sigma=dum$CVAR)
	}else{	
	    dum <- DrawLambdaPhi(RT,theta,sigma2,muI,SigmaI,ingroup)
	    if(Discrimination){	
			ab[,1] <- dum$phi  
			ab[,1] <- ab[,1]/(prod(ab[,1])^(1/K))
		}else{
			ab[,1] <- rep(1,K)  
	    }
		ab[,2] <- dum$lambda
	}
    
    MAB[ii,1:K,1:2] <- ab
    sigma2 <- SampleS_LNRT(RT,theta,ab[,2],ab[,1],ingroup)
	if(WL == 1){	    
		Msigma2[ii,1:K] <- 1/sqrt(sigma2) 	
	}else{
		Msigma2[ii,1:K] <- sigma2 	
	}
  
    X <- matrix(1,N,1)
    muP <- SampleB_LNRT(theta,X,SigmaP,muP0,SigmaP0)
    MmuP[ii,] <- muP$B
    muP <- muP$pred
   
    SS <- crossprod(theta - muP) + SigmaP0
    SigmaP <- rwishart_LNRT(1 + N, chol2inv(chol(SS)))$IW
    MSP[ii,,] <- SigmaP
    
	X <- matrix(1,K,1)
	if(WL == 1){
		ab1 <- cbind(1/sqrt(sigma2),ab[,2])
	}else{
		ab1 <- ab	
	}			
    muI2 <- SampleB_LNRT(ab1,X,SigmaI,muI0,SigmaI0)
    MmuI[ii,1] <- muI2$B[1]
    MmuI[ii,2] <- muI2$B[2]
    muI[,1] <- muI2$pred[,1]
    muI[,2] <- muI2$pred[,2]
    
    SS <- crossprod(ab1 - muI) + SigmaI0
    SigmaI <- rwishart_LNRT(2 + K, chol2inv(chol(SS)))$IW
    MSI[ii,,] <- SigmaI
    
    if(ii > 1000){  
		EAPphi <- (ab[,1] + (iis-1)*EAPphi)/iis
		EAPlambda <- (ab[,2] + (iis-1)*EAPlambda)/iis
		EAPtheta <- (theta + (iis-1)*EAPtheta)/iis
		EAPsigma2 <- (sigma2 + (iis-1)*EAPsigma2)/iis

		dum <- personfitLN(RT=RT,theta=theta,phi=ab[,1],lambda=ab[,2],sigma2=sigma2)
		lZP <- lZP + dum$lZP
		lZPT <- lZPT + dum$lZPT
		CF <- ifelse(dum$lZP < .05,1,0) #significance level = .05
		EAPCP <- (CF + (iis-1)*EAPCP)/iis

		dum <- itemfitLN(RT=RT,theta=theta,phi=ab[,1],lambda=ab[,2],sigma2=sigma2)
		lZI <- lZI + dum$lZI

		dum <- residualLN(RT=RT,theta=theta,phi=ab[,1],lambda=ab[,2],sigma2=sigma2,
			EAPtheta=EAPtheta,EAPlambda=EAPlambda,EAPphi=EAPphi,EAPsigma2=EAPsigma2)
		EAPresid <- EAPresid + dum$presid
		EAPKS <- (dum$KS[1:K,1] + (iis-1)*EAPKS)/iis

		iis <- iis + 1        
     }

	if(ii %% 100 == 0) cat("Iteration ",ii," ","\n")
	flush.console()
    }
  
  MT <- MT/XG
  MT2 <- sqrt(MT2/XG - MT^2 ) 
  if(ii > 1000){ 
   lZP <- lZP/(XG-1000)
   lZPT <- lZPT/(XG-1000)
   lZI <- lZI/(XG-1000) 	
   EAPresid <- EAPresid/(XG-1000) 	
  }
  
if(XG > 1000){
  return(list(Mtheta=MT,MTSD =MT2,MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,MSI = MSI,
	lZP =lZP,lZPT=lZPT,Msigma2=Msigma2,theta=theta,sigma2=sigma2,lZI=lZI,EAPresid=EAPresid,EAPKS=EAPKS,RT=RT,EAPCP=EAPCP))
}else{
  return(list(Mtheta=MT,MTSD =MT2,MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,MSI = MSI,
	Msigma2=Msigma2,theta=theta,sigma2=sigma2,RT=RT))
}
}




DrawLambdaPhi <- function(RT,theta,sigma2,muI,SigmaI,ingroup){
  
  library(MASS) 
  K <- ncol(RT)
  N <- nrow(RT)
  invSigmaI <- solve(SigmaI)
  H <- matrix(c(-theta,rep(1,N)),ncol=2,nrow=N)*ingroup 
  varest <- solve(kronecker(diag(1/sigma2[1:K]),(t(H)%*%H)) +	
                    kronecker(diag(1,K),invSigmaI))
  meanest <- t((t(H)%*%RT)/(t(matrix(sigma2,nrow=K,ncol=2))) + 
                 matrix(t(muI[1,]%*%invSigmaI),ncol=K,nrow=2)) 
  meanest <- apply((matrix(rep(meanest,K),ncol=2*K)%*%varest) * 
                     t(kronecker(diag(K),c(1,1))),2,sum)
  lambdaphi <- mvrnorm(1,mu=meanest,Sigma=varest)	
  lambdaphi <- matrix(lambdaphi,ncol=2,nrow=K,byrow=RT)	
  
  set <- which(lambdaphi[,1] < .30)
  if(length(set) >= 1){ 	
	  lambdaphi[set,1] <- .30
  }

  return(list(phi=lambdaphi[,1],lambda = lambdaphi[,2]))
}


SampleS_LNRT <- function(RT,zeta,lambda,phi,ingroup){
  K <- ncol(RT)
  N <- nrow(RT) 	
  Nr <- sum(ingroup)
  ss0 <- 10
  ingroup <- matrix(ingroup,ncol=K,nrow=N)
  Z <- (RT + t(matrix(phi,K,N))*matrix(zeta,N,K) - t(matrix(lambda,K,N)))*ingroup
  sigma2 <- (apply(Z*Z,2,sum) + ss0)/rchisq(K,Nr)
  return(sigma2)
}


DrawZeta_LNRT <- function(RT,phi,lambda,sigma2,mu,sigmaz){
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- matrix(lambda,nrow=K,ncol=N) - t(RT) 
  X <- matrix(phi,K,1)
  sigma2inv <- diag(1/sigma2[1:K]) 
  vartheta <- (1/((t(phi)%*%sigma2inv)%*%phi + 1/sigmaz))[1,1]
  meantheta <- matrix(((t(phi)%*%sigma2inv)%*%Z + t(mu/sigmaz))*vartheta,ncol=1,nrow=N)
  zeta <- matrix(rnorm(N,mean=meantheta,sd=sqrt(vartheta)),ncol=1,nrow=N)
  
  return(zeta)
}

SampleB_LNRT <- function(Y,X,Sigma,B0,V0){
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x%solve(crossprod(X)))+solve(V0))
  Btilde <- Bvar%*%(solve(Sigma%x%diag(m))%*%matrix(crossprod(X,Y),ncol=1)+solve(V0)%*%matrix(B0,ncol=1))
  B <- Btilde + chol(Bvar)%*% matrix(rnorm(length(Btilde)),ncol=1)
  pred <- X %*% matrix(B,ncol=k)
  return(list(B=B,pred=pred))
}

rwishart_LNRT <- function (nu, V) 
{
  m = nrow(V)
  df = (nu + nu - m + 1) - (nu - m + 1):nu
  if (m > 1) {
    RT = diag(sqrt(rchisq(c(rep(1, m)), df)))
    RT[lower.tri(RT)] = rnorm((m * (m + 1)/2 - m))
  }
  else {
    RT = sqrt(rchisq(1, df))
  }
  U = chol(V)
  C = t(RT) %*% U
  CI = backsolve(C, diag(m))
  return(list(IW = crossprod(t(CI))))
}


SimulateRT <- function(RT,zeta,lambda,phi,sigma2,DT){
#assume MAR 
	N <- nrow(RT)
	K <- ncol(RT)
  	meanT <- t(matrix(lambda,K,N)) - t(matrix(phi,K,N))*(zeta%*% t(rep(1,K)))
	meanT <- matrix(meanT,ncol=1,nrow=N*K)
      sigmaL <- matrix(t(matrix(rep(sqrt(sigma2),N),nrow=K,ncol=N)),ncol=1,nrow=N*K)
	RT <- matrix(RT,ncol=1,nrow=N*K) 
	RT[which(DT==0)] <- rnorm(sum(DT==0),mean=meanT[which(DT==0)],sd=sigmaL[which(DT==0)])
	RT <- matrix(RT,ncol=K,nrow=N)
	
	return(RT)
}


DrawLambda <- function(RT,zeta,sigma2,mu,sigma){

#prior mu,sigma
  N		<- nrow(RT)
  K 		<- ncol(RT) 
  zetaT	<- matrix(zeta,ncol=K,nrow=N,byrow=F)
  RTs 	<- matrix(RT+zetaT,ncol=K,nrow=N)	

  XX		<- matrix(1,ncol=K,nrow=N)
  pvar	<- diag(t(XX)%*%XX)/sigma2 + 1/sigma[1,1]
  betahat   <- diag(t(XX)%*%RT)/sigma2
  mu		<- (betahat + mu/sigma[1,1])/pvar
  beta	<-  rnorm(K,mean=mu,sd=sqrt(1/pvar))
	
	return(beta)
}


Conditional <- function(kk,Mu,Sigma,Z){
  
  K <- ncol(Z)
  N <- nrow(Z)

  if (kk == 1){
    C <- matrix(Z[,2:K] - Mu[,2:K],ncol=(K-1))
    CMEAN <- Mu[,1] + Sigma[1,2:K] %*% solve(Sigma[2:K,2:K]) %*% t(C)
    CSD <- Sigma[1,1] - Sigma[1,2:K] %*% solve(Sigma[2:K,2:K]) %*% Sigma[2:K,1]
  }
  if (kk > 1) {
    if(kk < K){
      C <- matrix(Z[1:N,1:(kk-1)] - Mu[,1:(kk-1)],ncol=(kk-1))
      CMu1 <- Mu[,kk:K] + t(matrix(t(Sigma[1:(kk-1),kk:K]),ncol=(kk-1)) %*% solve(Sigma[1:(kk-1),1:(kk-1)]) %*% t(C))
      CSigma <- Sigma[kk:K,kk:K] - matrix(t(Sigma[1:(kk-1),kk:K]),ncol=(kk-1)) %*% solve(Sigma[1:(kk-1),1:(kk-1)]) %*% matrix((Sigma[1:(kk-1),kk:K]),nrow=(kk-1))
      J <- ncol(CSigma)
      C <- matrix(Z[1:N,(kk+1):K] - CMu1[1:N,2:J],ncol=(J-1))
      CMEAN <- CMu1[,1] + CSigma[1,2:J] %*% solve(CSigma[2:J,2:J]) %*% t(C)
      CSD <- CSigma[1,1] - CSigma[1,2:J] %*% solve(CSigma[2:J,2:J]) %*% CSigma[2:J,1]
    }	
    if(kk ==K) {
      C <- matrix(Z[1:N,1:(K-1)] - Mu[,1:(K-1)],ncol=(K-1))
      CMEAN <- Mu[,K] + t(Sigma[1:(K-1),K]) %*% solve(Sigma[1:(K-1),1:(K-1)]) %*% t(C)
      CSD <- Sigma[K,K] - t(Sigma[1:(K-1),K]) %*% solve(Sigma[1:(K-1),1:(K-1)]) %*% Sigma[1:(K-1),K]
    }
  }
  return(list(CMU = CMEAN, CVAR=CSD))
}
