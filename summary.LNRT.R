
#' summaryRT
#' 
#' summary function for LNRT, and LNGAMMA
#' 
#' @param out     estimation result
#' @param data    simulated data
#' @export


summaryRT <- function(out,data){

#out : estimates
#data : simulated

if(missing(data)){
	simv <- FALSE
} 
else{
	simv <- TRUE
}
  N <- length(out$Mtheta)
  K <- ncol(out$MAB[,,1])
  XG <- length(out$MAB[,1,1])

  if(sum(out$MAB[,,1])==XG*K){
	WL <- 1 #no time discrimination
  }else{
	WL <- 0  
  }


if("Mnug" %in% names(out)){
	gammamodel <- TRUE 
}
else{
	gammamodel <- FALSE
}

N <- length(out$Mtheta)
K <- ncol(out$MAB[,,1])
XG <- length(out$MAB[,1,1])
bi <- round(.10*XG,0)  #burnin

##item parameter estimates

tdiscr <- apply(out$MAB[bi:XG,,1],2,mean)
tintens <- apply(out$MAB[bi:XG,,2],2,mean)

setdiscr <- round(sqrt(apply(out$MAB[bi:XG,,1],2,var)),3)
setintens <- round(sqrt(apply(out$MAB[bi:XG,,2],2,var)),3)

##item population parameter estimates

pdiscr <- round(apply(out$MmuI[bi:XG,],2,mean),3)
sepdiscr <- round(sqrt(apply(out$MmuI[bi:XG,],2,var)),3)

pSdiscr <- c(round(apply(out$MSI[bi:XG,,1],2,mean),3),round(apply(out$MSI[bi:XG,,2],2,mean),3))
sepSdiscr <- c(round(sqrt(apply(out$MSI[bi:XG,,1],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,,2],2,var)),3))

##person population parameter estimates

ppers <- round(mean(out$MmuP[bi:XG,]),3)
seppers <- round(sqrt(var(out$MmuP[bi:XG,])),3)

pSpers <- round(mean(out$MSP[bi:XG,1,1]),3)
sepSpers <- round(sqrt(var(out$MSP[bi:XG,1,1])),3)


##Shape parameter Gamma
# single nu : dim(Mnug) : XG*1
# multiple nu : dim(Mnug) : XG*K

if(gammamodel){
  if (ncol(out$Mnug)==1){
    estnug <- mean(out$Mnug[bi:XG])
    seestnug <- round(sqrt(var(out$Mnug[bi:XG])),3)
  }
  if (ncol(out$Mnug)>1){
    estnug <- round(apply(out$Mnug[bi:XG,],2,mean),3)
    seestnug <- round(sqrt(apply(out$Mnug[bi:XG,],2,var)),3)
  }
	
}
else{
  ##Measurement error parameter estimates
  
  estsigma2 <- round(apply(out$Msigma2[bi:XG,],2,mean),3)
  seestsigma2 <- round(sqrt(apply(out$Msigma2[bi:XG,],2,var)),3)
  
}

##Person Fit outcomes

if(gammamodel){
	cat("\n", "Gamma RT Modeling, 2013, J.P. Fox")
}
else{
	cat("\n", "Log-Normal RT Modeling, 2013, J.P. Fox")
}
	cat("\n", "Summary of results")


if(gammamodel){
  if(simv){
    if (ncol(out$Mnug)==1){
      cat("\n\n\t","Time Intensity parameter","\n")
      cat("\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
    }
    if (ncol(out$Mnug)>1){
      cat("\n\n\t","Time Intensity parameter","\t","Shape Parameter Gamma","\n")
      cat("\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\t","item", "\t","EAP","\t","SD","\n")
    }    
  }
  else{
    
    if (ncol(out$Mnug)==1){
      cat("\n\n\t","Time Intensity parameter","\n")
      cat("\t","item", "\t", "EAP", "\t", "SD","\n")
    }
    
    if (ncol(out$Mnug)>1){
      cat("\n\n\t","Time Intensity parameter","\t","Shape Parameter Gamma","\n")
      cat("\t","item", "\t", "EAP", "\t", "SD","\t\t","item", "\t","EAP","\t","SD","\n")      
    }    
    
  }
}

else{
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
}}

if(gammamodel){
  
  for(ii in 1:K){
    cat("\n\t", ii,"\t") 
    printF(tintens[ii], w=6,d=3)  # estimated bsp  
    cat("\t") 
    printF(setintens[ii], w=6, d=3)   
    if(simv){
      cat("\t")
      printF(data$bsp[ii],w=6, d=3)  # bsp true value
      #cat("\n\t", ii,"\t") 
      if(ncol(out$Mnug)>1){
        cat("\t")
        cat("\t", ii,"\t") 
        printF(estnug[ii],w=6, d=3)  # Shape Parameter Gamma EAP
        cat("\t")
        printF(seestnug[ii],w=6, d=3)        
      }
      
    }
    else{
      if(ncol(out$Mnug)==1){
        cat("\t")
      }
      
      if(ncol(out$Mnug)>1){
        cat("\t\t",ii,"\t")
        printF(estnug[ii],w=6, d=3)
        cat("\t")
        printF(seestnug[ii],w=6, d=3)        
      }      
    }
  }
}

else{
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

}


cat("\n\n\t", "Mean and Covariance matrix Items (phi,lambda)", "\n")
cat("\n\t", "--- Population Mean Item ---", "\n")
cat("\t","mu_phi", " ", "SD","\t","mu_lambda"," ","SD","\n\t")
	for(jj in c(1,2)){
		printF(pdiscr[jj], w=6,d=3)
		cat("\t")
		printF(sepdiscr[jj], w=6,d=3)
		cat("\t")
	}
cat("\n\n\t", "--- Covariance matrix Items ---", "\n")
cat("\t","phi", "\t", "SD", "\t", "Cov","\t","SD","\t", "lambda"," ","SD","\n\t")
	for(jj in c(1,2,4)){
		printF(pSdiscr[jj], w=6,d=3)
		cat("\t")
		printF(sepSdiscr[jj], w=6,d=3)
		cat("\t")
	}


cat("\n\n\t", "Mean and Covariance matrix Persons", "\n")
cat("\n\t", "--- Population Mean Person ---", "\n")
cat("\t","mu_P", "\t", "SD","\n\t")
	for(jj in c(1)){
		printF(ppers[jj], w=6,d=2)
		cat("\t")
		printF(seppers[jj], w=6,d=2)
		cat("\t")
	}
cat("\n\n\t", "--- Covariance matrix Person ---", "\n")
cat("\t","Sigma_P", "\t","SD","\n\t")
	for(jj in c(1)){
		printF(pSpers[jj], w=6,d=3)
		cat("\t")
		printF(sepSpers[jj], w=6,d=3)
		cat("\t")
	}

if( (gammamodel) && (ncol(out$Mnug)==1) ){
	cat("\n\n\t", "--- Shape Parameter Gamma ---", "\n")
	cat("\t","EAP", "\t", "SD","\n\t")
	printF(estnug, w=6,d=3)
	cat("\t")
	printF(seestnug[jj], w=6,d=3)
	cat("\t")
}


if("lZP" %in% names(out)){

cat("\n\n\n\t", "*** Person Fit Analysis ***", "\n")

#Percentage Aberrant lZP3
cat("\n\t", "Percentage Outliers Persons (5% level)","\n")
cat("\n\t", "lZ","\n")
cat("\t",round(100*sum(out$lZP < .05)/N,2),"%","\t\n")
cat("\t 95% Posterior Probability: ",round(100*sum(out$EAPCP > .95)/N,2),"%","\t\n")


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
}
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