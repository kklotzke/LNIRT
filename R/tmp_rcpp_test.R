library("devtools")
library("LNIRT")
library("microbenchmark")
library("rbenchmark")


find_rtools()

n1 <- 10
k1 <- 4

a1 <- runif(k1)
b1 <- runif(k1)
g1 <- runif(k1)
t1 <- runif(n1)
mu1 <- runif(k1)
muI1 <- matrix(runif(k1*2), ncol = 2)
m1 <- matrix(runif(4), ncol = 2)
m2 <- matrix(runif(1))
y1 <- matrix(runif(n1*k1), ncol=k1, nrow=n1)
y1 <- ifelse(y1 < 0.5, 0, 1)
d1 <- matrix(runif(n1*k1), ncol=k1, nrow=n1)
d1 <- ifelse(y1 < 0.5, 0, 1)
s1 <- matrix(runif(n1*k1), ncol=k1, nrow=n1)
s1 <- ifelse(s1 < 0.5, 0, 1)
d1 <- matrix(runif(n1*k1), ncol=k1, nrow=n1)
d1 <- ifelse(d1 < 0.5, 0, 1)
sig1 <- abs(rnorm(1))
rt1 <- matrix(rnorm(n1*k1), ncol=k1, nrow=n1)
sig21 <- runif(k1)
sigI1 <- diag(runif(2))
nu1 <- n1 + 2
nu2 <- n1 + 1
ingroup1 <- rep(1, n1)
ingroup2 <- sample(0:1, n1, replace = T)
theta1 <- matrix(rnorm(n1 * 2), ncol = 2) # 1: ability, 2: speed
muP1 <- matrix(rnorm(n1 * 2), ncol = 2) # Mean estimates for person parameters 
SigmaP1 <- diag(rnorm(2,1,1)) # Person covariance matrix


### Common functions ###

## Conditional

# kk == 1
outSimC11 <- Rcpp_Conditional(kk = 1, Mu = muP1, Sigma = SigmaP1, Z = theta1)
outSimC12 <- LNIRT:::Conditional(kk = 1, Mu = muP1, Sigma = SigmaP1, Z = theta1)
print(outSimC11)
print(outSimC12)

# kk == k
outSimC21 <- Rcpp_Conditional(kk = 2, Mu = muP1, Sigma = SigmaP1, Z = theta1)
outSimC22 <- LNIRT:::Conditional(kk = 2, Mu = muP1, Sigma = SigmaP1, Z = theta1)
print(outSimC21)
print(outSimC22)

# kk < k

## SimulateRT
set.seed(1)
outSimRT1 <- Rcpp_SimulateRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, sigma2 = sig21, DT = d1)
benchmark(Rcpp_SimulateRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, sigma2 = sig21, DT = d1))
print(outSimRT1)

set.seed(1)
outSimRT2 <- LNIRT:::SimulateRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, sigma2 = sig21, DT = d1)
benchmark(LNIRT:::SimulateRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, sigma2 = sig21, DT = d1))
print(outSimRT2)

all.equal(outSimRT1, outSimRT2)


## SimulateY
set.seed(1)
outSimY1 <- Rcpp_SimulateY(Y = y1, theta = t1, alpha0 = a1, beta0 = b1, guess0 = g1, D = d1)
benchmark(Rcpp_SimulateY(Y = y1, theta = t1, alpha0 = a1, beta0 = b1, guess0 = g1, D = d1))
print(outSimY1)

set.seed(1)
outSimY2 <- LNIRT:::SimulateY(Y = y1, theta = t1, alpha0 = a1, beta0 = b1, guess0 = g1, D = d1)
benchmark(LNIRT:::SimulateY(Y = y1, theta = t1, alpha0 = a1, beta0 = b1, guess0 = g1, D = d1))
print(outSimY2)

all.equal(outSimY1, outSimY2)



## DrawZeta
set.seed(1)
outZe1 <- Rcpp_DrawZeta(RT = rt1, phi = a1, lambda = b1, sigma2 = sig21, mu = t1, sigmaz = sig1)
benchmark(Rcpp_DrawZeta(RT = rt1, phi = a1, lambda = b1, sigma2 = sig21, mu = t1, sigmaz = sig1))
print(outZe1)

set.seed(1)
outZe2 <- LNIRT:::DrawZeta(RT = rt1, phi = a1, lambda = b1, sigma2 = sig21, mu = t1, sigmaz = sig1)
benchmark(LNIRT:::DrawZeta(RT = rt1, phi = a1, lambda = b1, sigma2 = sig21, mu = t1, sigmaz = sig1))
print(outZe2)

all.equal(outZe1, outZe2)


## rwishart
set.seed(1)
outR1 <- Rcpp_rwishart(nu = nu2, V = m2)
benchmark(Rcpp_rwishart(nu = nu1, V = m1))
print(outR1)

set.seed(1)
outR2 <- LNIRT:::rwishart(nu = nu2, V = m2)
benchmark(LNIRT:::rwishart(nu = nu1, V = m1))
print(outR2)

all.equal(outR1, outR2$IW)








### LNRT ###

## DrawLambdaPhi
set.seed(1)
outLP1 <- Rcpp_DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2)
benchmark(Rcpp_DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2))
print(outLP1)

set.seed(1)
outLP2 <- LNIRT:::DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2)
benchmark(LNIRT:::DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2))
print(outLP2)

outLP1 <- Rcpp_DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2)
lptmp <- outLP1$lambda
for (i in 1:1000)
{
  outLP1 <- Rcpp_DrawLambdaPhi_LNRT(RT = rt1, theta = t1, sigma2 = sig21, muI = muI1, sigmaI = sigI1, ingroup = ingroup2)
  lptmp <- rbind(lptmp, outLP1$lambda)
}

outLP3 <- mvrnorm(1000, outLP1$mu, outLP1$sigma)
colMeans(lptmp)
colMeans(outLP3)

## SampleS
set.seed(1)
outS_1 <- Rcpp_SampleS_LNRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, ingroup = ingroup2)
benchmark(Rcpp_SampleS_LNRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, ingroup = ingroup1))
print(outS_1)

set.seed(1)
outS_2 <- LNIRT:::SampleS_LNRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, ingroup = ingroup2)
benchmark(LNIRT:::SampleS2_LNRT(RT = rt1, zeta = t1, lambda = a1, phi = b1, ingroup = ingroup1))
print(outS_2)

all.equal(outS_1[, 1], outS_2)


## DrawLambda
set.seed(1)
outL1 <- Rcpp_DrawLambda_LNRT(RT = rt1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = sig1)
benchmark(Rcpp_DrawLambda_LNRT(RT = rt1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = sig1))
print(outL1)

set.seed(1)
outL2 <- LNIRT:::DrawLambda_LNRT(RT = rt1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = matrix(sig1,1,1))
benchmark(LNIRT:::DrawLambda_LNRT(RT = rt1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = matrix(sig1,1,1)))
print(outL2)

all.equal(outL1[,1], outL2)















### LNIRT ###

## DrawS
set.seed(1)
out1 <- Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1)
benchmark(Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1))

#set.seed(1)
#out2 <- Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1, alt = 0)
#benchmark(Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1, alt = 0))

set.seed(1)
out3 <- LNIRT:::DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1)
benchmark(LNIRT:::DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1))

all.equal(out1$S, out3)


## DrawZ
set.seed(1)
outZ1 <- Rcpp_DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = s1, D = d1, eta = matrix(0), PNO = 0)
benchmark(Rcpp_DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = s1, D = d1, eta = matrix(0), PNO = 0))

set.seed(1)
outZ2 <- LNIRT:::DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = s1, D = d1)
benchmark(LNIRT:::DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = s1, D = d1))


## DrawS + DrawZ
set.seed(1)
out1 <- Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1)
set.seed(1)
Rcpp_DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = out1$S, D = d1, eta = matrix(0), PNO = 0)
set.seed(1)
Rcpp_DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = out1$S, D = d1, eta = out1$eta, PNO = 1)
bench1 <- function()
{
  out1 <- Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1)
  Rcpp_DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = out1$S, D = d1, eta = out1$eta, PNO = 1)
}

bench2 <- function()
{
  out1 <- LNIRT:::DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1)
  LNIRT:::DrawZ_LNIRT(alpha0 = a1, beta0 = b1, theta0 = t1, S = out1, D = d1)
}
benchmark(bench1())
benchmark(bench2())


## DrawTheta
set.seed(1)
outT1 <- Rcpp_DrawTheta_LNIRT(alpha0 = a1, beta0 = b1, Z = outZ1, mu = t1, sigma = sig1)
benchmark(Rcpp_DrawTheta_LNIRT(alpha0 = a1, beta0 = b1, Z = outZ1, mu = t1, sigma = sig1))
print(outT1)

set.seed(1)
outT2 <- LNIRT:::DrawTheta_LNIRT(alpha0 = a1, beta0 = b1, Z = outZ1, mu = t1, sigma = sig1)
benchmark(LNIRT:::DrawTheta_LNIRT(alpha0 = a1, beta0 = b1, Z = outZ1, mu = t1, sigma = sig1))
print(outT2)

all.equal(outT1[,1], outT2)


## DrawBeta
set.seed(1)
outB1 <- Rcpp_DrawBeta_LNIRT(theta = t1, alpha = a1, Z = outZ1, mu = mu1, sigma = sig1)
benchmark(Rcpp_DrawBeta_LNIRT(theta = t1, alpha = a1, Z = outZ1, mu = mu1, sigma = sig1))
print(outB1)

set.seed(1)
outB2 <- LNIRT:::DrawBeta_LNIRT(theta = t1, alpha = a1, Z = outZ1, mu = mu1, sigma = matrix(sig1,1,1))
benchmark(LNIRT:::DrawBeta_LNIRT(theta = t1, alpha = a1, Z = outZ1, mu = mu1, sigma = matrix(sig1,1,1)))
print(outB2)

all.equal(outB1[, 1], outB2)


## DrawLambda
set.seed(1)
outL1 <- Rcpp_DrawLambda_LNIRT(RT = rt1, phi = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = sig1)
benchmark(Rcpp_DrawLambda_LNIRT(RT = rt1, phi = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = sig1))
print(outL1)

set.seed(1)
outL2 <- LNIRT:::DrawLambda_LNIRT(RT = rt1, phi = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = matrix(sig1,1,1))
benchmark(LNIRT:::DrawLambda_LNIRT(RT = rt1, phi = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = matrix(sig1,1,1)))
print(outL2)

all.equal(outL1, outL2$lambda)


## DrawAlpha
set.seed(1)
outA1 <- Rcpp_DrawAlpha_LNIRT(theta = t1, beta = b1, Z = outZ1, mu = mu1, sigma = sig1)
benchmark(Rcpp_DrawAlpha_LNIRT(theta = t1, beta = b1, Z = outZ1, mu = mu1, sigma = sig1))
print(outA1)

set.seed(1)
outA2 <- LNIRT:::DrawAlpha_LNIRT(theta = t1, beta = b1, Z = outZ1, mu = mu1, sigma = matrix(sig1,1,1))
benchmark(LNIRT:::DrawAlpha_LNIRT(theta = t1, beta = b1, Z = outZ1, mu = mu1, sigma = matrix(sig1,1,1)))
print(outA2)

all.equal(outA1[, 1], outA2)


## DrawPhi
set.seed(1)
outP1 <- Rcpp_DrawPhi_LNIRT(RT = rt1, lambda = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = sig1)
benchmark(Rcpp_DrawPhi_LNIRT(RT = rt1, lambda = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigmal = sig1))
print(outP1)

set.seed(1)
outP2 <- LNIRT:::DrawPhi_LNIRT(RT = rt1, lambda = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = matrix(sig1,1,1))
benchmark(LNIRT:::DrawPhi_LNIRT(RT = rt1, lambda = a1, zeta = t1, sigma2 = sig21, mu = mu1, sigma = matrix(sig1,1,1)))
print(outP2)

all.equal(outP1, outP2)


## SampleS2
set.seed(1)
outS21 <- Rcpp_SampleS2_LNIRT(RT = rt1, zeta = t1, lambda = a1, phi = b1)
benchmark(Rcpp_SampleS2_LNIRT(RT = rt1, zeta = t1, lambda = a1, phi = b1))
print(outS21)

set.seed(1)
outS22 <- LNIRT:::SampleS2_LNIRT(RT = rt1, zeta = t1, lambda = a1, phi = b1)
benchmark(LNIRT:::SampleS2_LNIRT(RT = rt1, zeta = t1, lambda = a1, phi = b1))
print(outS22)

all.equal(outS21[, 1], outS22)


microbenchmark("Rcpp" = Rcpp_DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1), 
               "R" = LNIRT:::DrawS_LNIRT(alpha0 = a1, beta0 = b1, guess0 = g1, theta0 = t1, Y = y1))









if(alt) {
  arma::mat a0(N, K);
  arma::mat b0(N, K);
  arma::mat g0(N, K);
  arma::mat t0(N, K);
  
  for (int i = 0; i < N; i++) {
    a0.row(i) = alpha0.t();
    b0.row(i) = beta0.t();
    g0.row(i) = guess0.t();
  }
  
  for (int i = 0; i < K; i++) {
    t0.col(i) = theta0;
  }
  
  NumericVector tmp = as<NumericVector>(wrap(a0 % t0 - b0));
  arma::mat eta(pnorm(tmp));
  eta.reshape(N, K);
  arma::mat probS = eta / (eta + g0 % (1 - eta));
  probS.reshape(N * K, 1);
  
  arma::vec randS(runif(N * K));
  arma::mat S(N * K, 1);
  S.fill(0);
  
  for(int i = 0; i < N * K; i++) {
    if(randS(i) > probS(i))
      S.row(i) = 1;
  }
  S.reshape(N, K);
  
  return(S % Y);
}