\name{SSbrook}
\alias{SSbrook}
\title{Selfstarting function for estimating hydraulic parameters in the 
Brooks-Corey model  }
\description{This selfstarting function estimates the water retention 
hydraulic parameters contained in the Brooks-Corey model. 
It has an initial attribute that creates initial estimates of 
the parameters Thr, Ths, alp, and scal 
}
\usage{
SSbrook(input, Thr, Ths, alp, scal)
}
\arguments{
  \item{input}{A numeric vector of natural logarithm of the suction potentials.}
  \item{Thr}{A numeric parameter representing the residual moisture content.}
  \item{Ths}{A numeric parameter representing the saturated moisture content.}
  \item{alp}{A numeric parameter representing the inverse of air-entry potential.}
  \item{scal}{A numeric parameter representing the pore-size distribution index.}
}
\details{The dataset should contain a column for logarithm of 
suction potential, taking care to place a very small numer for log (x) whenever x is zero. 
The function estimates the parameters in two parts: when product of alp and input is 
less than or eaual to unity and when this quantity is geater than unity
}
\value{A numeric vector of the same length as input. It is the value of 
the expression 
  \eqn{ Thr+(Ths-Thr)/((input*alp)^scal)  }.
It represents the fitted estimates of moisture contents using the Brooks-Corey model  }
\references{Brooks, R.H., Corey, A.T., 1964
  \emph{ Hydraulic properties of porous medium}
  \bold{ Hydrology Paper, Nr. 3}
   Colorado State University: fort Collins  }
\author{Omuto, Christian Thine}
\note{If the model returns negative values for the residual moisture content, 
a zero value should be dopted for this parameter. One may need to
 adjust the nls controls incase of convergence problems.It is also important 
to note that some datasets do perform poorly with this function  }
\seealso{SSvan, SSgard, SScamp}
\examples{
##---- Obtaine the data----
x <-c(0,log(10),log(31),log(100),log(200),log(500),log(2500),log(15000))# get log (x)
y <-c(0.362,0.353,0.294,0.230,0.206,0.182,0.148,0.136)
Brookdata <-data.frame(x,y)# make data frame
Brook.nls <-nls(y~SSbrook(x, Thr, Ths, alp, scal), Brookdata)# model for whole dataset
summary(Brook.nls)}
\keyword{ models}
