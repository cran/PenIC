\name{EM_fit}
\alias{EM_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for fitting generalized odds-rate model with specified link function (i.e., alpha value) under interval-censored data
}
\description{
Fits the generalized odds-rate model based on penalized B-splines to interval censored data via an EM algorithm.
}
\usage{
EM_fit(g0,b0,d1,d2,d3,Li,Ri,Z,nsub,alpha,qn,order,t.seq,tol=1e-5,itmax=500,lamu=1e5)
}
\arguments{
  \item{g0}{initial estimate of the spline coefficients; should be of length qn+order+1.}
  \item{b0}{initial estimate of regression coefficients; should be of length dim(Z)[2].}
  \item{d1}{vector indicating whether an observation is left-censored (1) or not (0).}
  \item{d2}{vector indicating whether an observation is interval-censored (1) or not (0).}
  \item{d3}{vector indicating whether an observation is right-censored (1) or not (0).}
  \item{Li}{the left endpoint of the observed interval; if an observation is left-censored, its corresponding entry should be 0. }
  \item{Ri}{the right endpoint of the observed interval; if an observation is right-censored, its corresponding entry should be Inf. }
  \item{Z}{design matrix of predictor variables (in columns); should be specified without an intercept term.}
  \item{nsub}{size of observed dataset.}
  \item{alpha}{parameter of link function; alpha=0 for the PH model and alpha=1 for the PO model.}
  \item{qn}{the number of interior knots to be used; should not exceed square root of sample size.}
  \item{order}{the order of the basis functions; order=3 for cubic spline.}
  \item{tol}{the convergence criterion of the EM algorithm.}
  \item{t.seq}{an increasing sequence of points at which the cumulative baseline hazard function is evaluated.}
  \item{itmax}{maximum iterations of EM procedure.}
  \item{lamu}{upper limit of smoothing parameter.}
}
\details{
The above function fits the generalized odds-rate model (with specified value of alpha) to interval censored data via an EM algorithm using penalized monotone B-splines.
}
\value{
  \item{b}{estimates of the regression coefficients.}
  \item{g}{estimates of the spline coefficients.}
  \item{se}{the standard deviation of b.}
  \item{base}{estimated cumulative baseline hazard function evaluated at the points t.seq.}
  \item{lambda}{final value of smooth parameter.}
  \item{flag}{the indicator whether the procedure converged; 0 if converged.}
}

\references{
Lu, M., Liu, Y., Li, C. and Sun, J. (2019+). An efficient penalized estimation approach for
a semi-parametric linear transformation model with interval-censored data. arXiv:1912.11703.
}
\examples{
set.seed(1)
case  <- 2
nsub  <- 35

# Generate interval-censored data under PH model

dat <- dataPA(nsub,case,alpha=0)
rp  <- c(mean(dat$d1),mean(dat$d2),mean(dat$d3))
rp

# [1] 0.2571429 0.3428571 0.4000000

t.seq <- seq(0.01,4,0.01)

# number of interior knots to be used
qn    <- ceiling(nsub^(1/3))-2
order <- 3
d1    <- dat$d1
d2    <- dat$d2
d3    <- dat$d3
Ri    <- dat$Ri
Li    <- dat$Li
Z     <- dat$Z
p     <- ncol(Z)
b0    <- rep(0,p)
g0    <- sort(runif(qn+order+1,-1,1))

# Fit data under PH model

fit <- EM_fit(g0,b0,d1,d2,d3,Li,Ri,Z,nsub,alpha=0,qn,order,t.seq,tol=1e-2,itmax=100,lamu=1e5)
cbind(fit$b,fit$se)


#           [,1]      [,2]
#[1,] -1.0655212 0.5021835
#[2,]  0.7649178 0.3185045

}
