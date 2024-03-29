% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LNRT_Q.R
\name{LNRTQ}
\alias{LNRTQ}
\title{Log-normal response time modelling with variable person speed (intercept, trend, quadratic)}
\usage{
LNRTQ(RT, X, data, XG = 1000, burnin = 10)
}
\arguments{
\item{RT}{a Person-x-Item matrix of log-response times (time spent on solving an item).}

\item{X}{explanatory (time) variables for random person speed (default: (1:N.items - 1)/N.items).}

\item{data}{either a list or a simLNIRTQ object containing the response time matrix.
If a simLNIRTQ object is provided, in the summary the simulated time parameters are shown alongside of the estimates.
If the RT variable cannot be found in the list, or if no data object is given, then the RT variable is taken
from the environment from which LNRTQ is called.}

\item{XG}{the number of MCMC iterations to perform (default: 1000).}

\item{burnin}{the percentage of MCMC iterations to discard as burn-in period (default: 10).}
}
\value{
an object of class LNRTQ.
}
\description{
Log-normal response time modelling with variable person speed (intercept, trend, quadratic)
}
