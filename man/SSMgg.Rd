\name{SSMgg}
\alias{SSMgg}
\title{ A new function based on the
 biexponential model to estimate soil hydraulic paremeters of water retention  }
\description{
  This selfstarting function estimates the water retention hydraulic 
parameters contained in the new biexponential function. It has an initial 
attribute that creates initial estimates of the parameters Ths1, alp1, Ths2, and alp1  }
\usage{
SSMgg(input, Ths1, alp1, Ths2, alp2)
}
\arguments{
  \item{input}{A numeric vector of the suction potentials  }
  \item{Ths1}{ A numeric parameter representing the saturated moisture 
content in the structural pore-space (macropores)  }
  \item{alp1}{ A numeric parameter representing the inverse of air-entry potential in the structural component  }
  \item{Ths2}{A numeric parameter representing the saturated moisture
 content in the textural pore-space (micropores)  }
  \item{alp2}{A numeric parameter representing the inverse of air-entry potential in the textural component  }
}
\value{
  numeric vector of the same length as input. It is the value of the expression
  \eqn{ Ths1*exp(-exp(alp1)*input) + Ths2*exp(-exp(alp2)*input) }. 
It represents the fitted estimates of moisture using the new biexponential model  }
\references{Groenevelt, P.H. and Grant, C.D., 2004
  \emph{ A new model for the soil-water retention curve that solves the problem of residual water contents  } 
  \bold{ European J.ournal of Soil Science, 55, 479-485  }
}
\author{Omuto, Christian Thine}
\note{ Note that the parameters representing the inverse 
of air-entry potentials are log-transformed. These coeficients need to 
be back-transformed using exponential function  }
\seealso{SSvan, SSbrook, SScamp}
\examples{
##---- Load isric data----
data(isric)
nbiexp.nls <-nls(y~SSMgg(x, Ths1, alp1, Ths2, alp2), data=isric)
plot(nbiexp.nls, Country~resid(.), abline=0) # to observe effect of geographic locations
}
\keyword{ models}
