## Fox : Code LNIRT Update May 2015
## Include marginal person-fit tests personfitR()
## Use parameterization a(theta-b) and phi(lambda-zeta)


LNIRT <- function(RT,Y,XG,guess,par1,residual,WL,td,alpha,beta){
  ## Main Programm function to call, uses other functions below 
  ## Inputs:
  ## Y = response matrix of dim(N=persons,K=items)
  ## RT = log-response time matrix (time spent on solving an item) of dim(N=persons,K=items) 
  ## XG = number of XG iterations for the MCMC algorithm
  ## guess: optional variable to indicate if guessing parameters should be included in the IRT model 
  ## td: optional variable to indicate if time-discrimination should be set to one, if missing td is included
  ## alpha :: pre-defined discrimination parameters
  ## beta ::  pre-defined difficulty parameters

 ## ident = 1: Identification : fix mean item difficulty(intensity) and product item (time) discrimination responses and response times 
 ## ident = 2: Identification : fix mean ability and speed and product item discrimination responses and response times 
 ##
 ## ident <- 1 (default) 	
  ident <- 2 #(to investigate person fit using latent scores)

  if (missing(guess)) {
    PNO <- 0 
  }else{
 	PNO <- 1	#guessing included
  }

  if (missing(par1)) {
    par1 <- 0 
  }else{
 	par1 <- 1	#use different parameterization
  }

  if (missing(residual)) {
    residual <- 0 
  }else{
 	residual <- 1	#residual analysis included
  }

  if (missing(WL)) {
	WL <- 0 
  }else{
 	WL <- 1	#time discrimination = 1/sqrt(error variance)
  }

  if(missing(td)){
	td <- 0 
  }else{
	td <- 1	#time discrimination fixed to one
  }

  if(missing(alpha)){
	discr <- 0 #discrimination estimated
  }else{
	discr <- 1 #discrimination known
  }

  if(missing(beta)){
	diffc <- 0 #difficulties estimated
  }else{
	diffc <- 1 #difficulties known
  }
 
  N <- nrow(Y) 
  K <- ncol(Y) 	

## population theta (ability - speed)
  theta <- matrix(rnorm(N*2),ncol=2)
  muP <- matrix(0,N,2)
  SigmaP <-  diag(2)
  muP0 <- matrix(0,1,2)
  SigmaP0 <-  diag(2)

##population item (ability - speed)
  ab <- matrix(rnorm(K*4),ncol=4)
  ab[,c(1,3)] <- 1
  muI <- t(matrix(rep(c(1,0),2*K),ncol=K))
  muI0 <- muI[1,]
  SigmaI <- diag(4)
  SigmaI0 <- diag(4)/10
  for(ii in 1:4) {
	SigmaI0[ii,] <- SigmaI0[ii,]*rep(c(.5,3),2)
  }
  ifelse(PNO>0,guess0 <- rep(0.2,K),guess0 <- rep(0,K)) 

##storage
  MT <- MT2 <- array(0,dim=c(N,2))
  MAB <- array(0,dim=c(XG,K,4))
  Mguess <- matrix(0,XG,K)
  MmuP <- array(0,dim=c(XG,2))
  MmuI <- array(0,dim=c(XG,4))
  MSP <- array(0,dim=c(XG,2,2))
  MSI <- array(0,dim=c(XG,4,4))
  Msigma2 <- matrix(0,ncol=K,nrow=XG)
  sigma2 <- rep(1,K)	

  lZP <- 0
  lZPT <- 0
  lZI <- 0
  lZPAT <- matrix(0,ncol=1,nrow=N)		
  lZIA <- matrix(0,ncol=1,nrow=K)		
  lZPA <- matrix(0,ncol=1,nrow=N)		

  EAPl0 <- matrix(0,ncol=K,nrow=N)		
  PFl <- matrix(0,ncol=1,nrow=N)		
  PFlp <- matrix(0,ncol=1,nrow=N)		
  IFl <- matrix(0,ncol=1,nrow=K)		
  IFlp <- matrix(0,ncol=1,nrow=K)		

  EAPresid <- matrix(0,ncol=K,nrow=N)
  EAPresidA <- matrix(0,ncol=K,nrow=N)
  EAPKS <- matrix(0,ncol=1,nrow=K)
  EAPKSA <- matrix(0,ncol=1,nrow=K)

  EAPphi <- matrix(0,ncol=1,nrow=K)
  EAPlambda <- matrix(0,ncol=1,nrow=K)
  EAPtheta <- matrix(0,ncol=2,nrow=N)
  EAPsigma2	<- matrix(0,ncol=1,nrow=K)

	if(XG > 1000){
 		MPF <- matrix(0,ncol=N,nrow=XG-1000) 
  		MPFb <- matrix(0,ncol=N,nrow=XG-1000) 
 		MPFp <- matrix(0,ncol=N,nrow=XG-1000) 
 		MPFbp <- matrix(0,ncol=N,nrow=XG-1000) 

  		MPFTb <- matrix(0,ncol=N,nrow=XG-1000) 
 		MPFTbp <- matrix(0,ncol=N,nrow=XG-1000) 
	}
	EAPbeta <- rep(0,K)
	EAPalpha <- rep(0,K)
	EAPmub <- 0
	EAPsigmab <- 1
	EAPSigmaP <- 1
	EAPmuI <- 0
	EAPSigmaI <- 1
	iis <- 1

  D <- matrix(1,ncol=1,nrow=N*K) 
  D[which(is.na(Y))] <- 0 
  D <- matrix(D,nrow=N,ncol=K)

  DT 	<- matrix(1,ncol=1,nrow=N*K) 
  DT[which(is.na(RT))] <- 0 
  DT <- matrix(DT,nrow=N,ncol=K)

  EAPCP1  <- matrix(0,ncol=1,nrow=N)	
  EAPCP2  <- matrix(0,ncol=1,nrow=N)	
  EAPCP3  <- matrix(0,ncol=1,nrow=N)	

  ## Start MCMC algorithm
  
  for (ii in 1:XG){
	if(sum(DT==0) > 0){
		if(WL == 1){
			RT <- SimulateRT(RT=RT,zeta=theta[,2],lambda=ab[,4],phi=rep(1,K),sigma2=sigma2,DT=DT)
		}else{
			RT <- SimulateRT(RT=RT,zeta=theta[,2],lambda=ab[,4],phi=ab[,3],sigma2=sigma2,DT=DT)
		}
	}
	#ability test
	if(PNO >0)  {
		if(sum(D==0) > 0){
			Y <- SimulateY(Y=Y,theta=theta[,1],alpha0=ab[,1],beta0=ab[,2],guess0=guess0,D=D)
		}
		if(par1==1){
			SR <- DrawS(alpha0=ab[,1],beta0=ab[,2]*ab[,1],guess0=guess0,theta=theta[,1],Y=Y)
			ZR <- DrawZ(alpha0=ab[,1],beta0=ab[,2]*ab[,1],theta0=theta[,1],S=SR,D=D)
		}else{
			SR <- DrawS(alpha0=ab[,1],beta0=ab[,2],guess0=guess0,theta=theta[,1],Y=Y)
			ZR <- DrawZ(alpha0=ab[,1],beta0=ab[,2],theta0=theta[,1],S=SR,D=D)		
		}
      }
	else{	
		if(par1==1){
			ZR <- DrawZ(alpha0=ab[,1],beta0=ab[,2]*ab[,1],theta0=theta[,1],S=Y,D=D)
		}else{
			ZR <- DrawZ(alpha0=ab[,1],beta0=ab[,2],theta0=theta[,1],S=Y,D=D)		
		}
	}

    dum <- Conditional(1,muP,SigmaP,theta)
	if(par1==1){
 		theta[,1] <- DrawTheta(alpha0=ab[,1],beta0=ab[,2]*ab[,1],Z=ZR,mu=dum$CMU,sigma=dum$CVAR)
	}else{
		theta[,1] <- DrawTheta(alpha0=ab[,1],beta0=ab[,2],Z=ZR,mu=dum$CMU,sigma=dum$CVAR)
	}	
	if(ident==2){#rescale for identification
		theta[,1] <- theta[,1] - mean(theta[,1])
	}	

    dum <- Conditional(2,muP,SigmaP,theta) 
	if(par1==1){
		theta[,2] <- DrawZeta(RT=RT,phi=ab[,3],lambda=ab[,3]*ab[,4],sigma2=sigma2,mu=as.vector(dum$CMU),sigmaz=as.vector(dum$CVAR)) ## speed 
	}else{
		theta[,2] <- DrawZeta(RT=RT,phi=ab[,3],lambda=ab[,4],sigma2=sigma2,mu=as.vector(dum$CMU),sigmaz=as.vector(dum$CVAR)) ## speed 
	}

	if(ident==2){#rescale for identification
		theta[,2] <- theta[,2] - mean(theta[,2])
	}	

    MT[1:N,1:2] <- MT[1:N,1:2] + theta
    MT2[1:N,1:2] <- MT2[1:N,1:2] + theta^2 
      
	if(PNO >0)  {
		guess0 <- DrawC(S=SR,Y=Y)
	}
    	Mguess[ii,] <- guess0

	if(WL == 1){
		ab1 <- cbind(ab[,1],ab[,2],1/sqrt(sigma2),ab[,4])
	}else{
		ab1 <- ab	
	}			

	if(discr==0){
	    	dum <- Conditional(kk=1,Mu=muI,Sigma=SigmaI,Z=ab1)#discrimination
		if(par1==1){
			ab[,1] <- abs(DrawAlpha(theta=theta[,1],beta=ab[,2],Z=ZR,mu=dum$CMU,sigma=dum$CVAR))
		}else{
			ab[,1] <- abs(DrawPhi(RT=ZR,lambda=-ab[,2],zeta=-theta[,1],sigma2=rep(1,K),mu=dum$CMU,sigmal=dum$CVAR))
		}
	}else{
		ab[,1] <- alpha
	}
	ab[,1] <- ab[,1]/(prod(ab[,1])^(1/K))
    
	if(WL == 1){
		ab1 <- cbind(ab[,1],ab[,2],1/sqrt(sigma2),ab[,4])
	}else{
		ab1 <- ab	
	}			

	if(diffc==0){
		dum <- Conditional(kk=2,Mu=muI,Sigma=SigmaI,Z=ab1)#difficulty
		if(par1==1){
			ab[,2] <- DrawBeta(theta=theta[,1],alpha=ab[,1],Z=ZR,mu=dum$CMU,sigma=dum$CVAR)
		}else{
			ab[,2] <- -DrawLambda(ZR,-ab[,1],theta[,1],rep(1,K),dum$CMU,dum$CVAR)$lambda
		}
	}else{
		ab[,2] <- beta
	}	
	if(par1==1){
			#ab[,2] <- ab[,2]/ab[,1]	
		}
	if(ident==1){#rescale for identification
		ab[,2] <- ab[,2] - mean(ab[,2])
	}	
	
	if((WL == 1) | (td == 1) ){
		#no time discrimination, 1/(sqrt(error variance)) = discrimination on MVN prior 
		ab[,3] <- 1	
	}else{	
	      dum <- Conditional(3,muI,SigmaI,ab)#time discrimination
		#if(par1==1){
		#    ab[,3] <- abs(DrawPhi(RT,ab[,3]*ab[,4],theta[,2],sigma2,dum$CMU,dum$CVAR))
		#}else{
		    ab[,3] <- abs(DrawPhi(RT,ab[,4],theta[,2],sigma2,dum$CMU,dum$CVAR))
		#}	
		ab[,3] <- ab[,3]/(prod(ab[,3])^(1/K))
	}

	if(WL == 1){
		ab1 <- cbind(ab[,1],ab[,2],1/sqrt(sigma2),ab[,4])
	}else{
		ab1 <- ab	
	}			
    dum <- Conditional(kk=4,Mu=muI,Sigma=SigmaI,Z=ab1)#time intensity
    ab[,4] <- DrawLambda(RT,ab[,3],theta[,2],sigma2,dum$CMU,dum$CVAR)$lambda
	#if(par1==1){
	#	ab[,4] <- ab[,4]/ab[,3]	
	#}
	if(ident==1){#rescale for identification
	    ab[,4] <- ab[,4] - mean(ab[,4])
	}	

    MAB[ii,1:K,1:4] <- ab
    sigma2 <- SampleS2(RT=RT,zeta=theta[,2],lambda=ab[,4],phi=ab[,3])
	if(WL == 1){	    
		Msigma2[ii,1:K] <- 1/sqrt(sigma2) 	
	}else{
		Msigma2[ii,1:K] <- sigma2 	
	}
       
    X <- matrix(1,N,1)
    muP <- SampleB_LNIRT(theta,X,SigmaP,muP0,SigmaP0)
    MmuP[ii,] <- muP$B
    muP <- muP$pred

    SS <- crossprod(theta - muP) + SigmaP0
    SigmaP <- rwishart_LNIRT(2 + N, chol2inv(chol(SS)))$IW
    MSP[ii,,] <- SigmaP

    X <- matrix(1,K,1)
	if(WL == 1){
		ab1 <- cbind(ab[,1],ab[,2],1/sqrt(sigma2),ab[,4])
	}else{
		ab1 <- ab	
	}			
    muI2 <- SampleB_LNIRT(Y=ab1,X=X,Sigma=SigmaI,B0=muI0,V0=SigmaI0)
    if(ident == 2){	
		MmuI[ii,c(1,2,3,4)] <- muI2$B[c(1,2,3,4)]
		muI[,c(1,2,3,4)] <- muI2$pred[,c(1,2,3,4)]
    }else{
		MmuI[ii,c(1,3)] <- muI2$B[c(1,3)]
		muI[,c(1,3)] <- muI2$pred[,c(1,3)]
    }	
    muI1 <- matrix(muI,ncol=4,nrow=K,byrow=FALSE)  
    SS <- crossprod(ab1 - muI1) + SigmaI0
    SigmaI <- rwishart_LNIRT(4 + K, chol2inv(chol(SS)))$IW
    MSI[ii,,] <- SigmaI
  
	if(ii > 1000){
		EAPmuI <- (muI[1,4] + (iis-1)*EAPmuI)/iis
		EAPsigma2 <- (sigma2 + (iis-1)*EAPsigma2)/iis
		EAPlambda <- (ab[,4] + (iis-1)*EAPlambda)/iis
		EAPphi <- (ab[,3] + (iis-1)*EAPphi)/iis
		EAPSigmaI <- (SigmaI[4,4] + (iis-1)*EAPSigmaI)/iis
	      EAPtheta <- (theta + (iis-1)*EAPtheta)/iis

		EAPmub <- (muI[1,2] + (iis-1)*EAPmub)/iis
		EAPsigmab <- (SigmaI[2,2] + (iis-1)*EAPsigmab)/iis
		if(par1 == 1){
			EAPbeta <- (ab[,1]*ab[,2] + (iis-1)*EAPbeta)/iis
		}else{
			EAPbeta <- (ab[,2] + (iis-1)*EAPbeta)/iis
		}
		EAPalpha <- (ab[,1] + (iis-1)*EAPalpha)/iis
		EAPSigmaP <- (SigmaP[1,1] + (iis-1)*EAPSigmaP)/iis


##############################

	if(residual){

## IRT Fit Evaluation
		if(par1==1){
			beta1 <- ab[,1]*ab[,2]
		}else{
			beta1 <- ab[,2] 
		}	
		if(PNO >0){
			dum <- residualA(Z=ZR,Y=SR,theta=theta[,1],alpha=ab[,1],beta=beta1,
						EAPtheta=EAPtheta[,1],EAPalpha=EAPalpha,EAPbeta=EAPbeta)
		}else{
			dum <- residualA(Z=ZR,Y=Y,theta=theta[,1],alpha=ab[,1],beta=beta1,
						EAPtheta=EAPtheta[,1],EAPalpha=EAPalpha,EAPbeta=EAPbeta)
		}
	      EAPKSA <- (dum$KS[1:K,1] + (iis-1)*EAPKSA)/iis
	      EAPresidA <- (dum$presidA + (iis-1)*EAPresidA)/iis

		#lZPAT <- lZPAT + dum$lZPAT
		lZPA <- lZPA + dum$lZPA
		#lZIA <- lZIA + dum$lZIA
		CF2 <- ifelse(dum$PFlp < .05,1,0) #significance level = .05
		EAPCP2 <- (CF2 + (iis-1)*EAPCP2)/iis


		EAPl0 <- ((iis-1)*EAPl0 + dum$l0)/iis
		PFl <- ((iis-1)*PFl + dum$PFl)/iis
		IFl <- ((iis-1)*IFl + dum$IFl)/iis
		PFlp <- ((iis-1)*PFlp + dum$PFlp)/iis
		IFlp <- ((iis-1)*IFlp + dum$IFlp)/iis

##############################

## Log-Normal Fit Evaluation
		if(par1 == 1){
			lambda1 <- ab[,3]*ab[,4]
		}else{
			lambda1 <- ab[,4]
		}				
	      dum <- personfitLN(RT=RT,theta=theta[,2],phi=ab[,3],lambda=lambda1,sigma2=sigma2)	# lZ statistic
	      lZP <- lZP + dum$lZP
	      lZPT <- lZPT + dum$lZPT
		CF1 <- ifelse(dum$lZP < .05,1,0) #significance level = .05
		EAPCP1 <- (CF1 + (iis-1)*EAPCP1)/iis #speed
		EAPCP3 <- (CF1*CF2 + (iis-1)*EAPCP3)/iis

		dum <- itemfitLN(RT=RT,theta=theta[,2],phi=ab[,3],lambda=lambda1,sigma2=sigma2)
	      lZI <- lZI + dum$lZI
	
		dum <- residualLN(RT=RT,theta=theta[,2],phi=ab[,3],lambda=lambda1,sigma2=sigma2,
			EAPtheta=EAPtheta[,2],EAPlambda=EAPlambda,EAPphi=EAPphi,EAPsigma2=EAPsigma2)
		EAPresid <- EAPresid + dum$presid
	      EAPKS <- (dum$KS[1:K,1] + (iis-1)*EAPKS)/iis
	  
        	iis <- iis + 1
	}
	
	}

	if(ii %% 100 == 0) cat("Iteration ",ii," ","\n")
	flush.console()

  }
	
  
  MT <- MT/XG
  MT2 <- sqrt(MT2/XG - MT^2 ) 
  
  if(ii > 1000){ 
	if(residual){
		lZP <- lZP/(XG-1000)
   		lZPT <- lZPT/(XG-1000)
		lZI <- lZI/(XG-1000) 	
		EAPresid <- EAPresid/(XG-1000) 	
		lZPAT <- lZPAT/(XG-1000)
		lZPA <- lZPA/(XG-1000)
		lZIA <- lZIA/(XG-1000) 	
  }	
}
  
 if(XG > 1000){
	if(residual){
  return(list(Mtheta=MT,MTSD =MT2, MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,MSI = MSI,Mguess=Mguess,Msigma2=Msigma2,
              lZP =lZP,lZPT=lZPT,lZPA=lZPA,lZI=lZI,EAPresid=EAPresid,EAPresidA=EAPresidA,EAPKS=EAPKS,EAPKSA,EAPKSA,
			PFl=PFl,PFlp=PFlp,IFl=IFl,IFlp=IFlp,EAPl0=EAPl0,RT=RT,Y=Y,EAPCP1=EAPCP1,EAPCP2=EAPCP2,EAPCP3=EAPCP3))
	}else{
	  return(list(Mtheta=MT,MTSD =MT2, MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,MSI = MSI,Mguess=Mguess,Msigma2=Msigma2,RT=RT,Y=Y))
	}
}else{
  return(list(Mtheta=MT,MTSD =MT2, MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,MSI = MSI,Mguess=Mguess,Msigma2=Msigma2,RT=RT,Y=Y))
}
}


DrawS <- function(alpha0,beta0,guess0,theta0,Y){
  N  		<- nrow(Y)
  K		<- ncol(Y) 
  eta 	<- t(matrix(alpha0,ncol=N,nrow=K)) * matrix(theta0,ncol=K,nrow=N)-t(matrix(beta0,nrow=K,ncol=N))
  eta 	<- matrix(pnorm(eta),ncol=K,nrow=N)
  probS 	<- eta/(eta + t(matrix(guess0,nrow=K,ncol=N))*(matrix(1,ncol=K,nrow=N)-eta))	
  S		<- matrix(runif(N*K),ncol=K,nrow=N)
  S		<- matrix(ifelse(S > probS,0,1),ncol=K,nrow=N)
  S		<- S*Y
  return(S)
}


DrawZ   <- function(alpha0,beta0,theta0,S,D){
  N 		<- nrow(S)
  K 		<- ncol(S) 
  eta	<- t(matrix(alpha0, ncol = N, nrow= K)) * matrix(theta0,ncol=K,nrow=N) - 
					t(matrix(beta0, ncol = N, nrow = K))
  BB		<- matrix(pnorm(-eta),ncol = K, nrow = N)
  BB[which(BB < .00001)] <- .00001
  BB[which(BB > (1-.00001))] <- (1 - .00001)
  u		<- matrix(runif(N*K), ncol = K, nrow = N)
  tt		<- matrix( ( BB*(1-S) + (1-BB)*S )*u + BB*S, ncol = K, nrow = N)
  


  Z  <- matrix(0,ncol=1,nrow=N*K)
  tt <- matrix(tt,ncol=1,nrow=N*K) 
  eta <- matrix(eta,ncol = 1, nrow = N*K) 
  Z[which(D==1)] <-  qnorm(tt)[which(D==1)] + eta[which(D==1)]
  Z[which(D==0)] <-  rnorm(sum(D==0)) + eta[which(D==0)]
  Z <- matrix(Z,ncol=K,nrow=N)

  return(Z)
}


DrawTheta <- function(alpha0,beta0,Z, mu, sigma) {
  
  N			<- nrow(Z)
  K 			<- ncol(Z) 
  pvar		<- (sum(alpha0^2)+1/as.vector(sigma))
  thetahat	<- (((Z + t(matrix(beta0,ncol= N, nrow = K))) %*% matrix(alpha0,ncol = 1,nrow = K)))
  mu			<- (thetahat +as.vector(mu)/as.vector(sigma))/pvar
  theta		<- rnorm(N,mean=mu,sd=sqrt(1/pvar))
  
  return(theta)
}


DrawZeta <- function(RT,phi,lambda,sigma2,mu,sigmaz){
  
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


DrawC <- function(S,Y){
  
  N		<- nrow(Y)
  K		<- ncol(Y) 
  Q1 		<- 20 + apply((S==0)*(Y==1),2,sum) #answer unknown and guess correctly
  Q2		<- 80 + apply(S==0,2,sum) #answer unknown
  guess	<- rbeta(K,Q1,Q2)
  return(guess)
}

DrawBeta <- function(theta,alpha,Z,mu,sigma){

#prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z) 
  alphaT	<- matrix(-alpha,ncol=K,nrow=N,byrow=T)
  thetaT	<- matrix(theta,ncol=K,nrow=N,byrow=F)
  Z 		<- matrix(Z-alphaT*thetaT,ncol=K,nrow=N)	

  XX		<- alphaT
  pvar	<- diag(t(XX)%*%XX) + 1/sigma[1,1]
  betahat   <- diag(t(XX)%*%Z)
  mu		<- (betahat + mu/sigma[1,1])/pvar
  beta	<-  rnorm(K,mean=mu,sd=sqrt(1/pvar))
	
	return(beta)
}


DrawLambda <- function(RT,phi,zeta,sigma2,mu,sigmal){
  
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- RT + t(matrix(phi,K,N))*matrix(zeta,N,K)
  X <- matrix(1,N,1) 
  Sigma <- diag(sigma2) 
  Sigma0 <- diag(K)*as.vector(sigmal)
  lambda <- SampleB_LNIRT(Z,X,Sigma,as.vector(mu),Sigma0)$B
  return(list(lambda = lambda))
}

DrawAlpha <- function(theta,beta,Z,mu,sigma){

#prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z) 
  betaT	<- matrix(beta,ncol=K,nrow=N,byrow=T)
  XX		<- matrix(c(theta-betaT),nrow=N,ncol=K)
  pvar	<- diag(t(XX)%*%XX) + 1/sigma[1,1]
  alphahat  <- diag(t(XX)%*%Z)
  mu		<- (alphahat + mu/sigma[1,1])/pvar
  alpha	<-  rnorm(K,mean=mu,sd=sqrt(1/pvar))
	
	return(alpha)
}

DrawPhi <- function(RT,lambda,zeta,sigma2,mu,sigmal){
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- -RT + t(matrix(lambda,K,N))
  X <- matrix(zeta,N,1) 
  Sigma <- diag(sigma2) 
  Sigma0 <- diag(K)*as.vector(sigmal)
  phi <- SampleB_LNIRT(Z,X,Sigma,as.vector(mu),Sigma0)$B
  return(phi)
}


SampleS2 <- function(RT,zeta,lambda,phi){
  K <- ncol(RT)
  N <- nrow(RT)
  ss0 <- 10
  Z <- RT + t(matrix(phi,K,N))*matrix(zeta,N,K) - t(matrix(lambda,K,N))
  sigma2 <- (apply(Z*Z,2,sum) + ss0)/rchisq(K,N)
  return(sigma2)
}


SampleB_LNIRT <- function(Y,X,Sigma,B0,V0){
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x%solve(crossprod(X)))+solve(V0))
  Btilde <- Bvar%*%(solve(Sigma%x%diag(m))%*%matrix(crossprod(X,Y),ncol=1)+solve(V0)%*%matrix(B0,ncol=1))
  B <- Btilde + chol(Bvar)%*% matrix(rnorm(length(Btilde)),ncol=1)
  pred <- X %*% matrix(B,ncol=k)
  return(list(B=B,pred=pred))
}

rwishart_LNIRT <- function (nu, V) 
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

SimulateY <- function(Y,theta,alpha0,beta0,guess0,D){
#with guessing
	N <- nrow(Y)
	K <- ncol(Y)

	G <- matrix(0,ncol=K,nrow=N)

	for(kk in 1:K){
		G[,kk] <- rbinom(N,size=1,prob=guess0[kk])
	}
	Y[(D==0) & (G==1)] <- 1 #missing: guessed correctly
 
	par <- theta %*% matrix(alpha0,nrow=1,ncol=K)-t(matrix(beta0,nrow=K,ncol=N))
	probs <- matrix(pnorm(par),ncol=K,nrow=N)	
  	Yn  <-matrix(runif(N*K),nrow = N, ncol = K)
	Yn <- ifelse(Yn < probs,1,0)
	Y[(D==0) & (G==0)] <- Yn[(D==0) & (G==0)] #missing: response generated

	return(Y)
}
	
		
