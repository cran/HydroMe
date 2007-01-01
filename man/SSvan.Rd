\name{SSvan}
\alias{SSvan}
\title{A selfstarting function for estimating hydraulic parameters 
in the van Genuchten model }
\description{
  This selfstarting function estimates the water retention parameters 
contained in the van Genuchtens function. It has an initial attribute that 
creates initial estimates of the parameters Thr, Ths, alp, and scal  }
\usage{
SSvan(input, Thr, Ths, alp, scal)
}
\arguments{
  \item{input}{A numeric vector of natural logarithm of the suction potentials  }
  \item{Thr}{A numeric parameter representing the residual moisture content  }
  \item{Ths}{A numeric parameter representing the saturated moisture content  }
  \item{alp}{A numeric parameter representing the inverse of air-entry potential  }
  \item{scal}{A numeric parameter representing the pore-size distribution index  }
}
\details{
 The dataset should contain a column for logarithm of suction 
potential, taking care to place zero for log (x) whenever x is zero. 
The function estimates the parameters in two parts: when the product of alp and input is less than 
or eaual to unity and when this quantity is geater than unity  }
\value{
  a numeric vector of the same length as input. 
It is the value of the expression 
 \eqn{ Thr+(Ths-Thr)/(1+(input*exp(alp))^scal)^(1-(1/scal)) }.  
It represents the fitted estimates of moisture retained using the van Genuchtens model  }
\references{ Van Genuchten, M.T., 1980
  \emph{ A closed-form Equation for predicting the hydraulic conductivity of unsaturated soils  } 
  \bold{ Soil Science Society of America Journal 44, 892-898  }
}
\author{ Omuto, Christian Thine}
\note{This model relates to a case of van Genuchtens model where
 m is 1 less than the inverse of n. The parameter for inverse of air-entry potential
 has been transformed using natural logarithm. Thus, the coefficient for this 
parameter needs to be transformed back using exponential function. 
It is also important to note that some datasets do perform poorly with this function  }
\seealso{ SSgard, SSbrook, SScamp}
\examples{
data(isric)
van.nls <-nls(y~SSvan(x, Thr, Ths, alp, scal), isric)
coef(van.nls)## Obtain the parameter estimates
##       Thr        Ths        alp       scal 
## 0.1747737  0.5010839 -2.9837774  1.3029102 

}
\keyword{ models}
