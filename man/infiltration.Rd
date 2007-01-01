\name{infiltration}
\alias{infiltration}
\docType{data}
\title{Infiltration characteristics of Awach watershed in
Lake Victoria Basin in Western Kenya}
\description{
  Part of infiltration characteristics from a study of soil physical properties 
in the Lake Victoria basin in Western Kenya}
\usage{data(infiltration)}
\format{
  A data frame with 4536 observations on the following 5 variables.
  \describe{
    \item{\code{Sample}}{a numeric vector of number of rows in the data  }
    \item{\code{Obj}}{an ordered factor with levels of sampling points  }
    \item{\code{Time}}{a numeric vector of time of infiltration  }
    \item{\code{Rate}}{a numeric vector of infiltration capacity  }
    \item{\code{Cumrate}}{a numeric vector of cummulative infiltration capacities  }
  }
}
\source{Omuto, C.T. 2003 
	\emph{Rapid Mapping of saturated hydraulic conductivity 
	in tropical watersheds  }
	\bold{MSc Thesis  }
	 Nairobi: University of Nairobi 
}
\examples{
data(infiltration)
## summary(infiltration) ; plot(infiltration) ...
}
\keyword{datasets}
