\name{isric}
\alias{isric}
\docType{data}
\title{Water retention characteristics of the ISRIC world databse}
\description{Part of the water retention data from several countries in the world 
as contained in the ISRIC world soil database. It contains names of the countries and sampling points as well
}
\usage{data(isric)}
\format{
  A data frame with 320 observations on the following 6 variables.
  \describe{
    \item{\code{Sample}}{a factor with levels with codes for samples  } 
    \item{\code{Country}}{a factor with levels with codefor countries  }
    \item{\code{BD}}{a numeric vector of subsoil bulk densities  }
    \item{\code{x}}{a numeric vector of suction potential in cm  }
    \item{\code{lnx}}{a numeric vector of logarithm of suction potential  }
    \item{\code{y}}{a numeric vector of moisture content in vol/vol  }
  }
}
\details{The data for each point contain at least five entries for water retention
 charcateristics including zero and 15000 cm suction potentials  }
\source{www.isric.nl  }
\examples{
data(isric)
## maybe summary(isric) ; plot(isric)
}
\keyword{datasets}
