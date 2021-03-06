\name{TBclustr-package}
\alias{TBclustr-package}
\alias{TBclustr}
\docType{package}
\title{
  TBclustr: Estimate reproductive numbers for tuberculosis from TB cluster sizes
}
\description{
  This package implements a model of importation and transmission of tuberculosis (TB) in a host population described in xxx paper. Importation of new cases is modelled as a poisson process, and the within-country transmission is modelled as a branching process with a given distribution of the number of offspring. Two types of distributions for secondary cases are considered: a negative binomial distribution and a log-normal distribution. The model is fitted to input data using EasyABC importation rate and within-population transmission rate to data specified in the cluster size file. 
}
\details{
modelstats function defines the branching process model. 

}
\author{
Leon Danon, l.danon@exeter.ac.uk
Ellen Brooks-Pollock, ellen.brooks-pollock@bristol.ac.uk
Maintainer: Leon Danon <l.danon@exeter.ac.uk>
}
\references{
  Brooks-Pollock et. al. in preparation, 2018. 
}
\keyword{ TBclustr,  }
\seealso{
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
