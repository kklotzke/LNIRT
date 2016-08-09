## fit stat in LNIRT

N <- 1000
K <- 10
# Simulate data: 
data <- simLNIRT(N=500,K=20,rho=0.9,WL=T)
data <- simLNIRT(N=500,K=20,rho=0.8)

out <- LNRT(RT=data$RT,XG=2000,WL=T) 
out1 <- LNRT(RT=data$RT,XG=2000) 
summaryRT(out=out1,data=data)
summaryRT(out=out,data=data)

plot(out2$PFl,out2$PFlp) #ability
plot(out2$EAPl0,out2$lZP,col="red")


out1 <- LNIRT(RT=data$RT1,Y=data$Y,XG=1000,WL=T) 
summaryIRT(out=out1,data)

out <- LNIRT(RT=data$RT,Y=data$Y,XG=2000,residual=T,WL=T) 
summaryIRT(out=out,data)


## Y = response matrix of dim(N=persons,K=items)
## T = log-response time matrix (time spent on solving an item) of dim(N=persons,K=items) 
## XG = number of iterations for the MCMC algorithm
## guess: optional variable to indicate if guessing parameters should be included in the IRT model. Use any number you like, e.g., guess = 1

out2 <- LNIRT(RT=data$RT,Y=data$Y,XG=1000)
summaryIRT(out=out2)
summaryIRT(out=out2,data=data)

out2 <- LNIRT(RT=data$RT,Y=data$Y,XG=1000,WL=T)
summaryIRT(out=out2)
summaryIRT(out=out2,data=data)


out21 <- LNIRT(RT=data$RT,Y=data$Yg,XG=1000,guess=T)
summaryIRT(out=out21,data=data)

out22 <- LNIRT(RT=data$RT1,Y=data$Y1,XG=1000,par1=T)
summaryIRT(out=out22,data=data)

out23 <- LNIRT(RT=data$RT1,Y=data$Y1g,XG=1000,par1=T,guess=T)
summaryIRT(out=out23,data=data)

out24 <- LNIRT(RT=data$RT1,Y=data$Y1g,XG=2000,par1=T,guess=T,residual=T)
summaryIRT(out=out24,data=data)

out25 <- LNIRT(RT=data$RT,Y=data$Y,XG=2000,residual=T)
summaryIRT(out=out25,data=data)

out26 <- LNIRT(RT=data$RT,Y=data$Yg,XG=2000,residual=T,guess=T)
summaryIRT(out=out26,data=data)

out27 <- LNIRT(RT=data$RT1,Y=data$Y1,XG=2000,residual=T,par1=T)
summaryIRT(out=out27,data=data)







out3 <- LNIRT(RT=data$RT2,Y=data$Y2,XG=1000,par1=T)
summaryIRT(out=out2)
summaryIRT(out=out3,data=data)

out4 <- LNIRT(RT=data$RT,Y=data$Y1,XG=2000,guess=T,residual=T,par1=T)
summaryIRT(out=out2)
summaryIRT(out=out4,data=data)


aa <- apply(out3$MAB[,,2],2,mean)
bb <- apply(out3$MAB[,,4],2,mean)


data <- simLNIRT(N=1000,K=20,rho=.8)
out <- LNIRT(Y=data$Y,RT=data$T,XG=2000)

Y <- data$Y
K <- ncol(data$Y)
sumscore <- apply(data$Y,1,sum)
set <- which(sumscore > K-1)
data$ab[order(data$ab[,2]),2]
Y[set,c(5,18)] <- 0 

out <- LNIRT(Y=Y,RT=data$T,XG=2000)

plot(out$Mtheta[,1],data$theta[,1])
plot(out$Mtheta[,2],data$theta[,2])

apply(out$MAB[,,1],2,mean)
apply(out$MAB[,,2],2,mean)
apply(out$MAB[,,3],2,mean)
apply(out$MAB[,,4],2,mean)
mean(out$MSP[,2,1]/(sqrt(mean(out$MSP[,1,1]))*sqrt(mean(out$MSP[,2,2]))))

plot(out$MAB[,2,1])
apply(out$MAB[,,1],2,mean)


plot(density(out$MPFp[1:1000,1]),ylim=c(0,4))
for(ii in 2:100){
  lines(density(out$MPFp[1:1000,ii]))
}
for(ii in 2:100){
  lines(density(out$MPFbp[1:1000,ii]),col="red")
}
for(ii in 2:100){
  lines(density(out1$MPFbp[1:1000,ii]),col="blue")
}
for(ii in set){
  lines(density(out$MPFp[1:1000,ii]),col="yellow")
}


plot(apply(out$MPF,2,mean),apply(out$MPFp,2,mean))
points(apply(out1$MPF,2,mean),apply(out1$MPFp,2,mean),col="red")

##points(apply(out$MPF,2,mean)[set],apply(out$MPFp,2,mean)[set],col="red")


plot(out$lZP1,apply(out$MPFp,2,mean))
points(out$lZP1[set],apply(out$MPFp,2,mean)[set],col="red")

plot(density(out1$lZP1))
lines(density(apply(out1$MPFp,2,mean)))


aa <- rnorm(K)

sum(aa%*%t(aa))
sum(aa)*sum(aa)


sum(aa**2) - sum(aa%*%t(aa))
sum((aa- mean(aa))**2)
sum(aa**2 - 2*mean(aa)*aa + mean(aa)**2)
sum(aa**2) - 2*mean(aa)*sum(aa) + K*mean(aa)**2
sum(aa**2) - 2*mean(aa)*sum(aa) + sum(aa)**2/K
sum(aa**2) - 2*sum(aa)**2/K + sum(aa)**2/K
sum(aa**2) - sum(aa)**2/K










