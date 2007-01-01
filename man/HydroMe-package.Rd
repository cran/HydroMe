\name{HydroMe-package}
\alias{HydroMe-package}
\alias{HydroMe}
\docType{package}
\title{A package for the estimation of soil hydraulic 
parameters using experimental data
}
\description{
This package uses experimental data for either infiltration or water retention 
and estimates hydraulic parameters contained in common infiltration and 
water retention models in soil science. The models considered include Philip's, 
Horton, and Green-Ampt among the infiltration functions and the van Genuchten, 
Brooks-Corey, Campbell, and Gardner's models among the water retention functions. 
The parameters of these functions are obtained by analytical curve-fitting method  
}
\details{The functions are selfstarting but use nlme library for parameter estimation. 
They are coded as follows: SSphilip, SShort, and SSGampt, respectively for the Philip's, 
Horton's, and Green-Ampt's infiltration functions and SSvan, SSbrook, 
SScamp, and SSgard, respectively for the van Genuchten's, Brooks-Corey's, Campbell's,
Gardner's, and new Biexponential models for water retention functions. 
The van Genuchten model is one in which m is 1 less the inverse of n
}
\author{Omuto, Christian Thine

Maintainer: Omuto Christian Thine <thineomuto@yahoo.com>
}
\references{ 
Omuto, C.T. and Gumbe, L.O.,  2007
   \emph{ R codes for estimating water infiltration and retention parameters using mixed-effects  }
   \bold{ Computers and Geosciences -Submitted } 
}
\keyword{ package }