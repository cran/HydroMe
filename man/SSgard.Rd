\name{SSgard}
\alias{SSgard}
\title{A selfstarting function to estimate hydraulic parameters in the
Gardners retention model}
\description{This selfstarting function estimates the 
water retention parameters contained in the Gardners function. 
It has an initial attribute that creates initial estimates of 
the parameters Thr, Ths, alp, and scal  }
\usage{
SSgard(input, Thr, Ths, alp, scal)
}
\arguments{
  \item{input}{A numeric vector of the suction potential  }
  \item{Thr}{A numeric parameter representing the residual moisture content  }
  \item{Ths}{A numeric parameter representing the saturated moisture content  }
  \item{alp}{A numeric parameter representing the inverse of air-entry potential  }
  \item{scal}{A numeric parameter representing the pore-size distribution index  }
}
\details{The dataset should contain a column for logarithm of 
suction potential, taking care to place a very small numer for log (x) whenever x is zero. 
The function estimates the parameters in two parts: when the product of alp and input is
 less than or eaual to unity and when this quantity is geater than unity  }
\value{numeric vector of the same length as input. It is the value of 
the expression 
 \eqn{ Thr+(Ths-Thr)/(1+(alp*input^scal)) }. 
It represents the fitted estimates of moisture contents using the Gardners model  }
\references{
 Gardner, W.R., 1958
 \emph{ Some steady state solutions of the unsaturated moisture flow equation 
with application to evaporation from a water table  }
 \bold{ Soil Science, 85, 228-232  }
}
\author{ Omuto, Christian Thine}
\note{If the model returns negative values for residual moisture 
content, a zero value should be dopted for that residual moisture content. 
Also, if the model fails to converge, one may need to adjust the controls. 
It is also important to note that some datasets do perform poorly with this function  }
\seealso{ SSvan, SSbrook, SScamp}
\examples{
##-- Load pfdata--
data(pfdata)
gard.nlis <-nlsList(y~SSgard(x, Thr, Ths, alp, scal)|Level, data=pfdata)
cor(fitted(gard.nlis),pfdata$y, use="complete")^2 # To compute r^2
}
\keyword{ models}
