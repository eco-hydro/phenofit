% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_param.R
\name{get_param}
\alias{get_param}
\alias{get_param.fFITs}
\title{Get parameters from curve fitting result}
\usage{
get_param(fits)

get_param.fFITs(fFITs)
}
\arguments{
\item{fits}{Multiple methods curve fitting results by \code{curvefits} result.}

\item{fFITs}{\code{fFITs} object returned by \code{\link[=curvefit]{curvefit()}}.}
}
\description{
Get parameters from curve fitting result
}
\examples{
library(phenofit)
# simulate vegetation time-series
fFUN = doubleLog.Beck
par  = c(
    mn  = 0.1,
    mx  = 0.7,
    sos = 50,
    rsp = 0.1,
    eos = 250,
    rau = 0.1)
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- fFUN(par, t)
methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
fFITs <- curvefit(y, t, tout, methods)
# multiple years
fits  <- list(`2001` = fFITs, `2002` = fFITs)

l_param   <- get_param(fits)
d_GOF     <- get_GOF(fits)
d_fitting <- get_fitting(fits)
l_pheno   <- get_pheno(fits, "AG", IsPlot=TRUE)
}