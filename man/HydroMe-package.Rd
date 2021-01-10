\name{HydroMe2-package}
\alias{HydroMe2-package}
\alias{HydroMe2}
\docType{package}
\title{
R codes for estimating water retention and infiltration model parameters using experimental data
}
\description{
This package is version 2 of HydroMe v.1 package. It estimates the parameters in infiltration and water retention models by curve-fitting method. The models considered are those that are commonly used in soil science. It has new models for water retention characteristic curve and debugging of errors in HydroMe  v.1
}
\details{
\tabular{ll}{
Package: \tab HydroMe2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-29\cr
License: \tab GPL\cr
}
The package contains hydraulic functions whose parameters are estimated from experimental data. The functions are supposed to be specified, the input data, and associated parameters declared in a standard regression modelling fashion. Where starting variables are needed, the package has a function known as Dstart for determining the starting variables
}
\author{
Christian Thine Omuto

Maintainer: Christian Thine Omuto <thineomuto@yahoo.com>

}
\references{
Omuto CT and Gumbe LO. 2009. Estimating water infiltration and retention characteristics using a computer program in R. Computers and Geosciences 35: 579-585
}

\keyword{ package }
\examples{
data(infilt)
require("minpack.lm")
gamp.ns <- nlsLM(Rate ~ SSgampt(Time,ks,A), infilt)
summary(gamp.ns)
data(isric)
isric1 <- isric[1:32,]
require("nlme")
omuto <- nlsList(y~SSomuto(x,ths1,alp1,ths2,alp2)|Sample, isric1)
omuto.nlme <- nlme(omuto)## for fitting mixed-effects models
summary(omuto.nlme)
}
