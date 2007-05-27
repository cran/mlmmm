\name{mlmmmbd.em}
\alias{mlmmmbd.em}
\title{
ML estimation under multivariate linear mixed models with block-diagonal covariance matrix and missing values
}
\description{
Implementation of em.pan() that restricts the covariance matrix
for the random effects to be block-diagonal. This function
is identical to pan() in every way except that psi is now 
characterized by a set of r matrices of dimension q x q.
}

\usage{
mlmmmbd.em(y, subj, pred, xcol, zcol, start, maxits=100, eps=0.01)
}
\arguments{
\item{y}{
See description for mlmmm.em().
}

\item{subj}{
See description for mlmmm.em().
}

\item{pred}{
See description for mlmmm.em().
}

\item{xcol}{
See description for mlmmm.em().
}

\item{zcol}{
See description for mlmmm.em().
}

\item{start}{
Same as for em.pan() except that the starting value for psi
have new dimensions: (q x q x r) 
}

\item{maxits}{
See description for mlmmm.em().
}

\item{eps}{
See description for mlmmm.em().

}}
\value{
A list with the same components as that from em.pan(), with a
minor difference: the dimension of "psi" is now (q x q x r).
}

\references{Schafer, J.L. and Yucel, R.M. (2002)
Computational strategies for multivariate linear mixed-effects
	     models with missing values, Journal of Computational and Graphical Statistics, 11, 421-442.

}

\keyword{models}

% Converted by Sd2Rd version 37351.
