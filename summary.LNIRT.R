
summaryIRT <- function(out,data){  
  
  if(missing(data)){
    simv <- FALSE
  } 
  else{
    simv <- TRUE
  }  
  if("Mnug" %in% names(out)){
    gammamodel <- TRUE 
  }
  else{
    gammamodel <- FALSE
  }
  
  N <- length(out$Mtheta[,1])
  K <- ncol(out$MAB[,,1])
  XG <- length(out$MAB[,1,1])
  bi <- round(.10*XG,0)  #burnin
 
  if(sum(out$MAB[,,3])==XG*K){
	WL <- 1 #no time discrimination
  }else{
	WL <- 0  
  }
  idiscr <- apply(out$MAB[bi:XG,,1],2,mean)   
  idiff <- apply(out$MAB[bi:XG,,2],2,mean) 
  tdiscr <- apply(out$MAB[bi:XG,,3],2,mean)
  tintens <- apply(out$MAB[bi:XG,,4],2,mean)

  seidiscr <- round(sqrt(apply(out$MAB[bi:XG,,1],2,var)),3) 
  seidiff <- round(sqrt(apply(out$MAB[bi:XG,,2],2,var)),3)
  setdiscr <- round(sqrt(apply(out$MAB[bi:XG,,3],2,var)),3)
  setintens <- round(sqrt(apply(out$MAB[bi:XG,,4],2,var)),3)  

  if(sum(out$Mguess)==0){
	Quess <- FALSE
  }else{
	Quess <- TRUE
      guess <- round(apply(out$Mguess,2,mean),3)
	seguess <- round(sqrt(apply(out$Mguess,2,var)),3)
  }	
  pdiscr2 <- round(apply(out$MmuI[bi:XG,],2,mean),3)    
  sepdiscr2 <- round(sqrt(apply(out$MmuI[bi:XG,],2,var)),3)
  
  pSdiscr2 <- matrix(c(round(apply(out$MSI[bi:XG,1,],2,mean),3),round(apply(out$MSI[bi:XG,2,],2,mean),3),
						round(apply(out$MSI[bi:XG,3,],2,mean),3),round(apply(out$MSI[bi:XG,4,],2,mean),3)),ncol=4,nrow=4,byrow=TRUE)
  sepSdiscr2 <- matrix(c(round(sqrt(apply(out$MSI[bi:XG,1,],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,2,],2,var)),3),
						round(sqrt(apply(out$MSI[bi:XG,3,],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,4,],2,var)),3)),ncol=4,nrow=4,byrow=TRUE)
  sds <- sqrt(diag(pSdiscr2))
  SigmaIcor <-  round(pSdiscr2/(sds%*%t(sds)),3)
  
  ppers2 <- round(apply(out$MmuP[bi:XG,],2,mean),3)
  seppers2 <- round(sqrt(apply(out$MmuP[bi:XG,],2,var)),3)
 
  pSpers2 <- matrix(c(round(mean(out$MSP[bi:XG,1,1]),3),round(mean(out$MSP[bi:XG,2,1]),3),round(mean(out$MSP[bi:XG,1,2]),3),
		round(mean(out$MSP[bi:XG,2,2]),3)),ncol=2,nrow=2)
  sepSpers2 <- matrix(c(round(sqrt(var(out$MSP[bi:XG,1,1])),3),round(sqrt(var(out$MSP[bi:XG,2,1])),3),
		round(sqrt(var(out$MSP[bi:XG,1,2])),3),round(sqrt(var(out$MSP[bi:XG,2,2])),3)),ncol=2,nrow=2)
  sds <- sqrt(diag(pSpers2))
  SigmaPcor <-  round(pSpers2/(sds%*%t(sds)),3)

  estsigma2 <- round(apply(out$Msigma2[bi:XG,],2,mean),3)
  seestsigma2 <- round(sqrt(apply(out$Msigma2[bi:XG,],2,var)),3)
 
  if(gammamodel){
    estnug <- mean(out$Mnug[bi:XG])
    seestnug <- round(sqrt(var(out$Mnug[bi:XG])),3)
  }
  if(gammamodel){
    cat("\n", "Gamma RT-IRT Modeling, 2013, J.P. Fox")
  }
  else{
    cat("\n", "Log-Normal RT-IRT Modeling, 2013, J.-P. Fox")
  }
  cat("\n", "Summary of results")
  
  if(simv){
    cat("\n\n\t", "Item Discrimination parameter", "\t", "Item Difficulty parameter","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t", "Sim","\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
  }
  else{
    cat("\n\n\t", "Item Discrimination parameter", "\t", "Item Difficulty parameter","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t\t","item", "\t", "EAP", "\t", "SD","\n")
  }
    for(ii in 1:K){
	cat("\t")
	if(simv){
		cat("\n\t", ii,"\t") 
		printF(idiscr[ii], w=6, d= 3) # EAP
		cat("\t")
		printF(seidiscr[ii],w=6, d=3) # SD
		cat("\t")
      	printF(data$ab[ii,1],w=6, d=3) # SIM
      	cat("\t",ii,"\t")
		printF(idiff[ii], w=6,d=3) 
		cat("\t") 
		printF(seidiff[ii], w=6, d=3)              	
		cat("\t")
		printF(data$ab[ii,2],w=6, d=3) 
	}else{ 
		cat("\n\t", ii,"\t") 
		printF(idiscr[ii], w=6, d= 3) # EAP
		cat("\t")
		printF(seidiscr[ii],w=6, d=3) # SD
	      cat("\t\t",ii,"\t")
    		printF(idiff[ii], w=6,d=3) 
		cat("\t") 
		printF(seidiff[ii], w=6, d=3)         
	 }
}  
  
if(simv){
	if(WL == 1){
		cat("\n\n\t","Time Discrimination (Measurement Error)"," ","Time Intensity","\n")
		cat("\t", "item", "\t", "EAP", "\t", "SD","\t", "Sim","\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
	}else{
		cat("\n\n\t","Time Discrimination","\t\t","Time Intensity","\t\t","Measurement Error Variance","\n")
		cat("\t", "item", "\t", "EAP", "\t", "SD","\t", "Sim","\t","item", "\t", "EAP", "\t", "SD","\t","Sim", "\t", "EAP", "\t", "SD", "\t", "Sim","\n")
	}
}
else{
	if(WL == 1){
		cat("\n\n\t","Time Discrimination (Measurement Error)","\t","Time Intensity","\n")
		cat("\t", "item", "\t", "EAP", "\t", "SD","\t\t\t","item", "\t", "EAP", "\t", "SD","\t","\n")
	}else{
		cat("\n\n\t","Time Discrimination","\t","Time Intensity","\t","Measurement Error Variance","\n")
		cat("\t","item", "\t", "EAP", "\t", "SD","\t","item", "\t", "EAP", "\t", "SD","\t", "EAP", "\t", "SD","\n")
	}	
}
for(ii in 1:K){  
	cat("\t")
	if(simv){
		cat("\n\t", ii,"\t") 
		if(WL == 1){ 
			printF(estsigma2[ii],w=6, d=3) # EAP
	  		cat("\t")
	  		printF(seestsigma2[ii],w=6, d=3) # SD     
			cat("\t")
	    		printF(1/sqrt(data$sigma2[ii]),w=6, d=3) #SIM      
		}else{
			printF(tdiscr[ii], w=6,d=3) # EAP
			cat("\t")
			printF(setdiscr[ii], w=6, d=3) # SD
			cat("\t")
			printF(data$ab[ii,3],w=6, d=3) #SIM   
		}
		cat("\t",ii,"\t")    
  		printF(tintens[ii], w=6,d=3) # EAP
  		cat("\t") 
  		printF(setintens[ii], w=6, d=3) # SD
    		cat("\t")
    		printF(data$ab[ii,4],w=6, d=3) # SIM      
		cat("\t")
		if(WL == 1){ 
  			cat(" ")
		}else{
			printF(estsigma2[ii],w=6, d=3) # EAP
  			cat("\t")
  			printF(seestsigma2[ii],w=6, d=3) # SD     
  			cat("\t")
	    		printF(data$sigma2[ii],w=6, d=3) #SIM      
		}
	}else{
		cat("\n\t", ii,"\t") 
		if(WL == 1){ 
			printF(estsigma2[ii],w=6, d=3) # EAP
	  		cat("\t")
	  		printF(seestsigma2[ii],w=6, d=3) # SD     
	  		cat("\t\t")
		}else{
			printF(tdiscr[ii], w=6,d=3) # EAP
			cat("\t")
			printF(setdiscr[ii], w=6, d=3) # SD
		}
		cat("\t",ii,"\t")    
  		printF(tintens[ii], w=6,d=3) # EAP
  		cat("\t") 
  		printF(setintens[ii], w=6, d=3) # SD
		cat("\t")
		if(WL == 1){ 
			cat("")
		}else{	
			printF(estsigma2[ii],w=6, d=3) # EAP
  			cat("\t")
  			printF(seestsigma2[ii],w=6, d=3) #SD     
		}
	}  
}



if(Quess){

if(simv){
  cat("\n\n\t","Guessing Parameter","\n")
  cat("\t", "item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
	for(ii in 1:K){
	  cat("\n\t", ii,"\t") 
	  # Guessing Parameter
	  printF(guess[ii], w=6,d=3) # EAP
	  cat("\t") 
	  printF(seguess[ii], w=6, d=3) # SD
	  cat("\t") 
	  printF(data$quess[ii], w=6, d=3) # true value
	}
}else{
  cat("\n\n\t","Guessing Parameter","\n")
  cat("\t", "item", "\t", "EAP", "\t", "SD","\n")
	for(ii in 1:K){
	  cat("\n\t", ii,"\t") 
	  # Guessing Parameter
	  printF(guess[ii], w=6,d=3) # EAP
	  cat("\t") 
	  printF(seguess[ii], w=6, d=3) # SD
	}
}
}

if(round(pdiscr2[2],3)==0 && round(pdiscr2[4],3)==0 && (WL == 0)){
 cat("\n\n\t", "Mean and Covariance matrix Items (mu_a,mu_phi)", "\n")
 cat("\n\t", "--- Population Mean Item ---", "\n")
 cat("\t","mu_a","\t", "SD","\t","mu_phi"," ","SD","\n")
 for(jj in c(1,3)){
    cat("\t")
    printF(pdiscr2[jj], w=6,d=3)
    cat("\t")
    printF(sepdiscr2[jj], w=6,d=3)
  }
}
else{
    cat("\n\n\t", "Mean and Covariance matrix Items (mu_a,mu_b,mu_phi,mu_lambda)", "\n")
    cat("\n\t", "--- Population Mean Item ---","\n")
    cat("\t","mu_a","\t", "SD","\t","mu_b","\t","SD","\t","mu_phi","\t", "SD","\t","mu_lambda","\t","SD","\n\t")
    for(jj in c(1,2,3,4)){
      printF(pdiscr2[jj], w=6,d=3)
      cat("\t")
      printF(sepdiscr2[jj], w=6,d=3)
      cat("\t")
    }    
}
if(WL==1){
  cat("\n\n\t", "--- Covariance matrix Items (a,b,error variance,lambda)---", "\n")
  cat("\t\t","SigmaI","\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
}else{
  cat("\n\n\t", "--- Covariance matrix Items (a,b,phi,lambda)---", "\n")
  cat("\t\t","SigmaI","\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
}
 for(ii in 1:4){
  cat("\t") 
  printF(pSdiscr2[ii,], w=6,d=3)
  cat("\t")
	if(max(sepSdiscr2)<1e02){
	  printF(sepSdiscr2[ii,], w=6,d=3)
	}else{
	 printF(sepSdiscr2[ii,], w=8,d=3)
	}
  cat("\t") 
  printF(SigmaIcor[ii,], w=6,d=3)
  cat("\t\n") 
 }  
 
 cat("\n\n\t", "Mean and Covariance matrix Persons (ability,speed)", "\n")
 cat("\n\t", "--- Population Mean Person (Ability - Speed)---", "\n")
 cat("\t\t","muP", "\t\t", "SD","\n")
 for(ii in 1:2){
  if(ii == 1) {
	cat("\t Ability ") 
  }else{	
	cat("\t Speed \t") 
  }
  printF(ppers2[ii], w=6,d=3)
  cat("\t") 
  printF(seppers2[ii], w=6,d=3)
  cat("\t\n") 
 }  
 cat("\n\t","SigmaP", "\t", "SD SigmaP","\t","SigmaP (Correlation)","\n")
 for(ii in 1:2){
  cat("\t") 
  printF(pSpers2[ii,], w=6,d=3)
  cat("\t") 
  printF(sepSpers2[ii,], w=6,d=3)
  cat("\t") 
  printF(SigmaPcor[ii,], w=6,d=3)
  cat("\t\n") 
 }  

if(gammamodel){
    cat("\n\n\t", "--- Shape Parameter Gamma ---", "\n")
    cat("\t","EAP", "\t", "SD","\n\t")
    printF(estnug, w=6,d=3)
    cat("\t")
    printF(seestnug[jj], w=6,d=3)
    cat("\t")
}

if("lZP" %in% names(out)){

cat("\n\n\n\t", "*** \t\t\t\t\t\t ***", "\n")
cat("\t", "*** Person Fit Analysis (Log-Normal Speed) ***", "\n")
cat("\t", "*** \t\t\t\t\t\t ***", "\n")


#Percentage Aberrant lZP (chi-square statstic)
cat("\n\t", "Percentage Outliers Persons (5% level)","\n")
cat("\n\t", "lZ","\n")
cat("\t",round(100*length(which(out$lZP < .05))/N,2),"%","\t\n")
cat("\t 95% Posterior Probability: ",round(100*sum(out$EAPCP1 > .95)/N,2),"%","\t\n")

# Percentage Persons (Estimation Group, Aberrant Group)

#	cat("\n\t", "Percentage Persons (Estimation Group, Aberrant Group) (5% level)","\n")
#	cat("\n\t", "Estimation Group", "\t", "Aberrant","\n")
#	cat("\t",100*round(apply(out$Mingroup,2,mean)[1],3),"\t\t\t",100*round(apply(out$Mingroup,2,mean)[2],3),"\t\n")


cat("\n\n\t", "*** Item Fit Analysis ***", "\n")

cat("\n\t", "Misfitting Items (5% level)","\n")
set <- which(out$lZI < .05)
if(length(set >= 1)){
	cat("\n\t", "lI","\n")
	cat("\t",set,"\t\n")
}else{
	cat("\t","No Misfitting Items","\t\n")
}

cat("\n\t", "*** Residual Analysis ***","\n")
set <- which(out$EAPresid > .95)
if(length(set >= 1)){
	cat("\n\t", "Percentage Extreme Residuals (.95 Posterior Probability)","\n")
	cat("\t",round(100*length(set)/N,4),"%","\t\n")

	## Identify Extremes 
	set <- which(out$EAPresid > .95)
	set <- cbind((set %% N),(floor(set/N)+1),exp(out$RT[set]))
	dd <- order(set[,1])
	colnames(set) <- c("Person","Item", " RT ")
	rownames(set) <- rep(c(""),nrow(set))

	cat("\n\t","Extreme Residuals","\n")
	cat("\t","Person"," Item","\t"," RT ","\n")
	for(jj in dd){
		cat("\t")
		printF(set[jj,1], w=6,d=0)
		cat("\t")
		printF(set[jj,2], w=6,d=0)
		cat("\t")
		printF(set[jj,3], w=8,d=4)
		cat("\n")
	}
	cat("\n\t", "Kolmogorov Smirnov Test (5% level)","\n")
	cat("\t",round(100*length(which(out$EAPKS < .05)))/K,"%","of items has non-lognormally distributed residuals","\t\n")
	set <- which(out$EAPKS < .05)
	if(length(set) > 0){ 
		cat("\t","Item","\t P-value ","\n")
		for(jj in 1:length(set)){
			cat("\t")
			printF(set[jj], w=6,d=0)
			cat("\t")
			printF(out$EAPKS[set[jj]], w=6,d=3)
			cat("\n")
		}
	}
}else{
	cat("\t","No Extreme Residuals","\t\n")
	cat("\n\t", "Kolmogorov Smirnov Test (5% level)","\n")
	cat("\t",round(100*length(which(out$EAPKS < .05)))/K,"%","\t\n")
	set <- which(out$EAPKS < .05)
	if(length(set) > 0){ 
		cat("\t","Item","\t P-value ","\n")
		for(jj in 1:length(set)){
			cat("\t")
			printF(set[jj], w=6,d=0)
			cat("\t")
			printF(out$EAPKS[set[jj]], w=6,d=3)
			cat("\n")
		}
	}
}

cat("\n\n\n\t", "*** \t\t\t\t\t\t ***", "\n")
cat("\t", "*** Person Fit Analysis (IRT Model For Ability) ***", "\n")
cat("\t", "*** \t\t\t\t\t\t ***", "\n")

#Percentage Aberrant lZP (chi-square statstic)
#Percentage Aberrant l0 (log-likelihood statstic)
cat("\n\t", "Percentage Outliers Persons (5% level)","\n")
cat("\n\t", "Log-likelihood Statistic","\n")
#cat("\t",round(100*length(which(out$lZPA < .05))/N,2),"%","\t\n")
cat("\t",round(100*length(which(out$PFlp < .05))/N,2),"%","\t\n")
cat("\t 95% Posterior Probability: ",round(100*sum(out$EAPCP2 > .95)/N,2),"%","\t\n")
cat("\t 95% Posterior Probability (Ability and Speed): ",round(100*sum(out$EAPCP3 > .95)/N,2),"%","\t\n")


cat("\n\n\t", "*** Item Fit Analysis ***", "\n")

cat("\n\t", "Misfitting Items (5% level)","\n")
#set <- which(out$lZIA < .05)
set <- which(out$IFlp < .05)
if(length(set >= 1)){
	cat("\n\t", "lI","\n")
	cat("\t",set,"\t\n")
}else{
	cat("\t","No Misfitting Items","\t\n")
}

cat("\n\t", "*** Residual Analysis ***","\n")
set <- which(out$EAPresidA > .95)
if(length(set >= 1)){
	cat("\n\t", "Percentage Extreme Residuals (.95 Posterior Probability)","\n")
	cat("\t",round(100*length(set)/N,4),"%","\t\n")

	## Identify Extremes 
	set <- which(out$EAPresidA > .95)
	set <- cbind((set %% N),(floor(set/N)+1),out$Y[set], out$Mtheta[(set %% N)])
	dd <- order(set[,1])
	colnames(set) <- c("Person","Item", " Y ", " EAP theta")
	rownames(set) <- rep(c(""),nrow(set))

	cat("\n\t","Extreme Residuals","\n")
	cat("\t","Person"," Item","\t"," Response ", " EAP Theta", "\n")
	for(jj in dd){
		cat("\t")
		printF(set[jj,1], w=6,d=0)
		cat("\t")
		printF(set[jj,2], w=6,d=0)
		cat("\t")
		printF(set[jj,3], w=6,d=0)
		cat("\t")
		printF(set[jj,4], w=8,d=4)
		cat("\n")
	}
	cat("\n\t", "Kolmogorov Smirnov Test (5% level)","\n")
	cat("\t",round(100*length(which(out$EAPKSA < .05)))/K,"%","of items has non-normally distributed latent residuals","\t\n")
	set <- which(out$EAPKSA < .05)
	if(length(set) > 0){ 
		cat("\t","Item","\t P-value ","\n")
		for(jj in 1:length(set)){
			cat("\t")
			printF(set[jj], w=6,d=0)
			cat("\t")
			printF(out$EAPKSA[set[jj]], w=6,d=3)
			cat("\n")
		}
	}
}else{
	cat("\t","No Extreme Residuals","\t\n")
	cat("\n\t", "Kolmogorov Smirnov Test (5% level)","\n")
	cat("\t",round(100*length(which(out$EAPKSA < .05)))/K,"%","\t\n")
	set <- which(out$EAPKSA < .05)
	if(length(set) > 0){ 
		cat("\t","Item","\t P-value ","\n")
		for(jj in 1:length(set)){
			cat("\t")
			printF(set[jj], w=6,d=0)
			cat("\t")
			printF(out$EAPKSA[set[jj]], w=6,d=3)
			cat("\n")
		}
	}
}
} ## close personfit report
cat("\n\n")
}



printF <- function(x, w, d, matrix = F){
if(matrix)
x <- as.matrix(x)
if(!is.numeric(x))
stop("x must be numeric in printF")
if(length(x) == 1) {
cat(formatF(x, w = w, d = d))
}
else if(!matrix) {
cat(formatF(as.vector(x), w = w, d = d))
}
else {
apply(as.matrix(x), 1, FUN = function(y, w, d)
{
cat(formatF(y, w = w, d = d))
}
, w = w, d = d)
}
invisible()
}

formatF <- function(x, w, d){

#
# Format x as a Fixed Point number with fixed-width w and d places to the 
#  right of the decimal pt. If the number is too wide to fit, fill
#  the field with '*'s.  If d=0, do not print the decimal point.
#  Pad left with blanks. Pad right with zeros.
#
#  x can be a vector.  All elements of the vector will be formatted
#   the same. x cannot be a matrix or data frame.
#
wholePart <- as.integer(x)
wholeStrings <- as.character(wholePart)
wholeStrings <- ifelse((wholeStrings == 0) & (x < 0), paste("-", 
wholeStrings, sep = ""), wholeStrings)
wholeWidths <- ifelse(d > 0, w - d - 1, w)
leftPad <- wholeWidths - nchar(wholeStrings)
decimalPart <- round(x - wholePart, digits = d)
decimalStrings <- as.character(decimalPart)
decimalStrings <- substring(decimalStrings, regexpr("\\.", 
decimalStrings) + 1)
rightPad <- ifelse(rep(d, length(x)) > 0, d - nchar(decimalStrings),
0)
for(i in seq(along = wholeStrings)) {
if(leftPad[i] >= 0) {
wholeStrings[i] <- paste(paste(rep(" ", leftPad[i]),
collapse = ""), wholeStrings[i], sep = "")
}
else {
wholeStrings[i] <- paste(rep("*", wholeWidths[i]),
collapse = "")
}
if(rightPad[i] > 0) {
decimalStrings[i] <- paste(decimalStrings[i], paste(
rep("0", rightPad[i]), collapse = ""), sep = 
"")
}
}
decimalPoint <- ifelse(d > 0, ".", "")
if(d > 0)
paste(wholeStrings, decimalPoint, decimalStrings, sep = "")
else paste(wholeStrings)
}