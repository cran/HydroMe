\name{Infilt}
\alias{Infilt}
\docType{data}
\title{Infiltration test data from the Upper Athi River basin, Eastern Kenya }
\description{
  
  This data contains infiltration capacities(in cm/min), Time (min), 
and land use-land cover types in the Upper Athi river basin in Eastern Kenya. 
The infiltration tests were carried out using standard infiltrometers and 
the data sampled over many plots in the whole basin.
}
\usage{data(Infilt)}
\format{
  A data frame with 1105 observations on the following 5 variables.
  \describe{
    \item{\code{Sample}}{a numeric vector of raws in the data  }
    \item{\code{PlotNo}}{an ordered factor of plocodes  }
    \item{\code{Erosion}}{a factor with levels \code{E0} \code{E1} \code{E2}  }
    \item{\code{Time}}{a numeric vector of the time of infiltration tests in minutes  }
    \item{\code{Rate}}{a numeric vector of infiltration capacities in cm/min  }
  }
}
\source{
  Omuto, C.T. 2006
    \emph{ Large-Area Soil Physical Degradation Assessment: The Use of GIS, Remote Sensing,and Infrared 
Spectroscopy in Arid and Semi-Arid Kenya}
    \bold{ PhD dissertation}
     Nairobi: University of Nairobi
}
\references{
 Omuto, C.T., Minasny, B., McBratney, A.B., Biamah, E.K., 2006
    \emph{ Nonlinear mixed effect modelling for improved estimation
 of water retention and infiltration parameters}
    \bold{ Journal of Hydrology, 330, 748-758}
}
\examples{
data(Infilt)
summary(Infilt) 
}
\keyword{datasets}