\name{dataPA}
\alias{dataPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Date generation function
}
\description{
Generate interval-censored data under generalized odds-rate model, with different combinations of right-censoring
rate and cumulative baseline hazard function.
}
\usage{
dataPA(N, case, alpha)
}
\arguments{
  \item{N}{size of dataset}
  \item{case}{data generation configuration; takes value in 1, 2 and 3. }
  \item{alpha}{parameter of link function; alpha=0 for the PH model and alpha=1 for the PO model.}
}
\details{
The above function generate interval-censored data from generalized odds-rate model, under different simulation configurations. For further details please see Lu et al. (2019+).
}
\value{
  \item{d1}{vector indicating whether an observation is left-censored (1) or not (0).}
  \item{d2}{vector indicating whether an observation is interval-censored (1) or not (0).}
  \item{d3}{vector indicating whether an observation is right-censored (1) or not (0).}
  \item{Li}{the left endpoint of the observed interval; if an observation is left-censored, its corresponding entry should be 0. }
  \item{Ri}{the right endpoint of the observed interval; if an observation is right-censored, its corresponding entry should be Inf. }
  \item{Z}{design matrix of predictor variables (in columns); should be specified without an intercept term.}
}

\references{
Lu, M., Liu, Y., Li, C. and Sun, J. (2019+). An efficient penalized estimation approach for a semi-parametric linear transformation model with interval-censored data. arXiv:1912.11703.
}
\examples{
case  <- 3
nsub  <- 100

# Generate interval-censored data under PH model

dat <- dataPA(nsub,case,alpha=0)
rp  <- c(mean(dat$d1),mean(dat$d2),mean(dat$d3))
rp

# [1] 0.63 0.22 0.15

}
