\name{SScamp}
\alias{SScamp}
\title{A selfstarting function for estimating hydraulic parameters in the
Campbell model}
\description{ This selfstarting model evaluates the Campbells 
water retention function and its parameters. It has an initial attribute 
that evaluates initial estimates of the parameters Ths, alp, and scal 
for a given set of data  }
\usage{
SScamp(input, Ths, alp, scal)
}
\arguments{
  \item{input}{numeric vector of natural logarithm of suction potential values at which to evaluate the model  }
  \item{Ths}{ a numeric parameter representing the saturated moisture content when suction potential is zero  }
  \item{alp}{a numeric parameter representing the inverse of air-entry potential  }
  \item{scal}{a numeric parameter representing the slope of the water retention curve. It is related to pore-size distribution index  }
}
\details{
 Campbell function is a three-parameter model that does not assume
 the residual moisture content. The function estimates the parameters 
in two parts: when the product of alp and input is less than or eaual to unity and when
 this quantity is geater than unity  }
\value{a numeric vector of the same length as input.It is the value of the Campbells 
expression 
 \eqn{ Ths*(alp*input)^(-scal) }.  
It represents the fitted estimates of moisture contents using the Campbell function  }
\author{Omuto, Christian Thine}
\seealso{ SSgard, SSbrook, SSvan}
\examples{
##---- use part of the Wret data ---- 
data(Wret)
nWret=nrow(Wret) # gets the number of rows need (there are 232 all together) 
testW13 <-Wret[(nWret-231):(nWret-208), ] # choose the first 3 sets 
Campbel.nls <-nls(y~SScamp(lnx, Ths, alp, scal), testW13) 
Campbel.nlis <-nlsList(y~SScamp(lnx, Ths, alp, scal)|Points, testW13) 
Campbel.nlme <-nlme(Campbel.nlis) # Estimate mixed-effects
}
\keyword{ models}