\name{pfdata}
\alias{pfdata}
\docType{data}
\title{ Water retention data from the Upper Athi river basin in Eastern Kenya}
\description{Part of water retention data obtained from a study of soil physical 
degradation in the Upper Athi river basin in Eastern Kenya
}
\usage{data(pfdata)}
\format{
  A data frame with 2196 observations on the following 6 variables.
  \describe{
    \item{\code{sample}}{a numeric vector of the number of rows in the data  }
    \item{\code{Level}}{an ordered factor with levels for sampling points  }
    \item{\code{Landuse}}{a factor with levels \code{Agroforestry} \code{Cropland} \code{Forest}
        \code{Grassland} \code{Shrubland} \code{Thicket}  }
    \item{\code{x}}{a numeric vector of suction potential in cm  }
    \item{\code{y}}{a numeric vector of moisture content in vol/vol  }
    \item{\code{lnx}}{a numeric vector of logarithm of suction potential  }
  }
}
\source{Omuto, C.T. 2006
	\emph{ Large-Area Soil Physical Degradation
 Assessment: The Use of GIS, Remote Sensing,and Infrared Spectroscopy
 in Arid and Semi-Arid Kenya}
	\bold{ PhD dissertation}
	Nairobi: University of Nairobi
}
\references{Omuto, C.T., Minasny, B., McBratney, A.B., Biamah, E.K., 2006
	\emph{ Nonlinear mixed effect modelling for improved estimation 
	of water retention and infiltration parameters }
	\bold{ Journal of Hydrology, 330, 748-758 }
}
\examples{
data(pfdata)
## maybe summary(pfdata) ; plot(pfdata) 
}
\keyword{datasets}