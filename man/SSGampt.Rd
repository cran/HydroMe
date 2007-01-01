\name{SSGampt}
\alias{SSGampt}
\title{A selfstarting function for estimating 
infiltration parameters in the Green-Ampt model}
\description{
  This selfstarting function estimates ks and (A) during infiltration.
 It has an initial attribute that creates initial estimates of the parameters ks and A  }
\usage{
SSGampt(input, ks, A)
}
\arguments{
  \item{input}{a numeric vector of values at which to 
evaluate the model. It should contain the cummulative infiltration capacities  }
  \item{ks}{ a numeric parameter representing the 
asymptote to the cummulative-infiltration axis (when input values are very large)  }
  \item{A}{a numeric parameter representing the integral 
measure of the potential of wetting front in increasing soil moisture during infiltration  }
}
\value{a numeric vector of the same length as input. It is the value of the
 expression 
 \eqn{ ks+(A/input) }. 
It represents 
the fitted estimates of infiltration rates  }
\references{Green, W.A., Ampt, G.A., 1911
 \emph{ Studies on soil physics: 1. the flow of air and water through soils}
 \bold{ Journal of Agricultural Science 4, 1-24  } 
}
\author{Omuto, Christian Thine }
\note{ It may be necessary to transform the instantaneous 
infiltration rates to avoid the function returning negative ks values }
\seealso{ SShort, SSphilip  }
\examples{
##---- Use infiltration data---- 
data(infiltration) 
Gampt.nls <-nls(log(Rate)~SSGampt(Cumrate, ks, A), data=infiltration) 
summary(Gampt.nls)# Call summary to desplay the results 

}
\keyword{ models}
