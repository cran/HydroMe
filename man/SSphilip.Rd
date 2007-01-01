\name{SSphilip}
\alias{SSphilip}
\title{A selfstarting function for estimating
 infiltration parameters in the Philips model  }
\description{This selfstarting function estimates the 
infiltration parameters contained in the Philips model. It has an initial attribute 
that creates initial estimates of the parameters fc, and S  }
\usage{
SSphilip(input, fc, S)
}
\arguments{
  \item{input}{ A numeric vector of Time  }
  \item{fc}{ A numeric parameter representing the final steady infiltration capacity  }
  \item{S}{ a numeric parameter representing the slope of the infiltration 
characteristic curve at the early tmies of the infiltration process. It is the estimate of sorptivity at the supply moisture content  }
}
\value{
  a numeric vector of the same length as input. It is the value of the 
expression 
 \eqn{ fc+(0.5*S)/sqrt(input) }. 
 It represents the fitted values by the Philips model }
\references{Philip, J.R., 1957
 \emph{ The theory of infiltration 4: Sorptivity and algebraic
 infiltration equations }
 \bold{ Soil Science 84, 257-264  }
}
\author{ Omuto, Christian Thine}
\seealso{SShort, SSGampt}
\examples{
##-- Use infiltration data 
data(infiltration)
 philip.nls <-nls(Rate~SSphilip(Time, fc, S), data=infiltration) 
plot(fitted(philip.nls)~infiltration$Rate) # Plot the fitted versus measured

}
\keyword{ models}
