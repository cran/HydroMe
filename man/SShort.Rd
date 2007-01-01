\name{SShort}
\alias{SShort}
\title{A selfstarting function for estimating the parameters of the Hortons infiltration
model  }
\description{This selfstarting function estimates the infiltration parameters 
contained in the Hortons model. It has an initial attribute that 
creates initial estimates of the parameters fc, f0, and k  }
\usage{
SShort(input, fc, f0, lrk)
}
\arguments{
  \item{input}{A numeric vector of Time  }
  \item{fc}{A numeric parameter representing the final steady infiltration capacity  }
  \item{f0}{A numeric parameter representing the initial infiltration capacity  }
  \item{lrk}{A numeric parameter representing the logarithm of the decay constant  }
}
\value{
  a numeric vector of the same length as input. It is the value 
of the expression 
  \eqn{ fc+(f0-fc)*exp(-input*lrk) }. 
It represents the fitted estimates of infiltration capacities using the Hortons model  }
\references{ Horton , R.E., 1940
  \emph{ An approach towards a physical interpretation of infiltration capacity  }
  \bold{ Soil Science Society of America Proceedings 5, 227-237  }
}
\author{Omuto, Christian Thine}
\note{ It may be necessary to transform the infiltration capacities to avoid 
negative values for the asymptotic final steady infiltration capacity  }
\seealso{SSphilip, SSGampt}
\examples{
##-- Use infiltration data 
data(infiltration) 
hort.nls <-nls(Rate~SShort(Time, fc, f0, lrk), data=infiltration) 
plot(hort.nls) # to view the assumption of constant variance of the residuals
}
\keyword{ models}
