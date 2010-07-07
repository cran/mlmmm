\name{mlmmm.em}
\alias{mlmmm.em}
\title{
ML estimation via EM-algorithm under multivariate linear mixed models with missing values
}
\description{
ML estimation via hybrid of EM and Fisher scoring algorithm under the multivariate linear mixed models
with missing values described by Schafer and Yucel (2002), Yucel (2007). This function
will typically be used to produce maximum likelihood estimation of the 
unknown parameters under the model

yi = Xi\%*\%beta + Zi\%*\%bi + ei ,    i=1,\dots,m,

where

yi    = (ni x r) matrix of incomplete multivariate
        data for subject or cluster i;

Xi    = (ni x p) matrix of covariates;

Zi    =	(ni x q) matrix of covariates;

beta  = (p x r) matrix of coefficients common to the
        population (fixed effects);

bi    = (q x r) matrix of coefficients specific to
        subject or cluster i (random effects); and

ei    = (ni x r) matrix of residual errors.


The matrix bi, when stacked into a single column, is assumed
to be normally distributed with mean zero and unstructured
covariance matrix psi, and the rows of ei are assumed to be
independently normal with mean zero and unstructured
covariance matrix sigma. Missing values may appear in yi in
any pattern.

In most applications of this model, the first columns of Xi
and Zi will be constant (one) and Zi will contain a subset of
the columns of Xi. 
}

\usage{
mlmmm.em(y, subj, pred, xcol, zcol, start, maxits=200, eps=0.0001)
}

\arguments{
\item{y}{
matrix of responses. This is simply the individual yi matrices
stacked upon one another. Each column of y corresponds to a
response variable. Each row of y corresponds to a single
subject-occasion, or to a single subject within a cluster.
Missing values (NA) may occur in any pattern.
}

\item{subj}{
vector of length nrow(y) giving the subject (or cluster)
indicators i for the rows of y. For example, suppose 
that y is in fact rbind(y1,y2,y3,y4) where nrow(y1)=2,
nrow(y2)=3, nrow(y3)=2, and nrow(y4)=7. Then subj should
be c(1,1,2,2,2,3,3,4,4,4,4,4,4,4).
}

\item{pred}{
matrix of covariates used to predict y. This should have the
same number of rows as y. The first column will typically be
constant (one), and the remaining columns correspond to other
variables appearing in Xi and Zi.
}

\item{xcol}{
vector of integers indicating which columns of pred will be
used in Xi. That is, pred[,xcol] is the Xi matrices (stacked
upon one another).
}

\item{zcol}{
vector of integers indicating which columns of pred will be
used in Zi. That is, pred[,zcol] is the Zi matrices (stacked
upon one another).
}

\item{start}{
optional list of quantities to specify the initial estimates
of the parameters for the EM. 
If "start" is omitted then mlmmm.em() chooses its own initial values.
}

\item{maxits}{
maximum number of cycles of EM to be performed. The algorithm runs 
to convergence or until "maxits" iterations, whichever comes first.
}

\item{eps}{
convergence criterion. The algorithm is considered to have
converged if the relative differences in all parameters from
one iteration to the next are less than eps--that is, if
all(abs(new-old)<eps*abs(old)).
}}

\value{
   A list containing  the following elements:
   \item{beta}{A matrix containing the final value of the estimate of the fixed effects. The first column corresponds to the estimates for the first column of y, the second column corresponds to the estimates of the second column  of y, and so on.
	       }

    \item{Sigma}{
A matrix containing the final value of the  estimate of the variance covariance matrix of the vectorized residual matrix term.
                }

\item{Psi}{
A matrix containing the final value of the estimate of the variance covariance matrix of the (vectorized) random-effects matrix. 
}

\item{eb}{
A matrix (of dimensions r*q by m) containing the emprical bayes estimates of the random-effects \code{b_i}.
} 

\item{varb}{
An array of dimensions r*q x r*r x m, containing the variance covariance matrix of the random-effects.
}

\item{xtwxinv}{
Variance-covariance matrix of the estimate of fixed estimates.
}

\item{converged}{
An indicator showing whether the algorithm converged or not.
}

\item{iter}{
Number of iterations to convergence.
}

\item{npatt}{
Number of distinct missingness patterns, not counting the ones missing all variables making the response matrix.
}

\item{pstfin}{
A matrix of dimensions npatt by r, indicating the number of rows with the underlying missingness pattern.
}

\item{iposn}{
A vector showing the row numbers of y, which belong to missingness patterns showed in pstfin.
}

\item{patt}{
A vector of n denoting the missingness patterns of the rows of y.
}

\item{rmat}{
A matrix showing the distinct missingness patterns, excluding the rows that are completely missing.
}

\item{logll}{
A vector of expected loglikelihood values at each iteration.
}

\item{logoll}{
A vector of observed loglikelihood values at each iteration.
}

\item{clock}{
How much time (in seconds) mlmmmm.em took to converge.
}}


\details{
The EM algorithm used in mlmmm.pan() is described in detail by
Schafer and Yucel (2002) and Yucel (2007). 
}

\references{
Schafer, J.L. and Yucel, R.M. (2002) Computational strategies for multivariate linear mixed-effects models with missing values. Journal of the Computational and Graphical Statistics, Volume 11, Number 2,
  437--457.

Yucel, R.M. (2007) R mlmmm package: Fitting multivariate linear mixed-effects models with missing values
}

\examples{\dontrun{
# For a detailed example, see the file "mlmmmmex.s" distributed
# with this function. Here is a simple example of how mlmmm.em()
# might be used to produce Ml estimates.
library(mlmmm)
data(adg)
y<-cbind(adg$y.1,adg$y.2)
colnames(y)=c("adg","initwt")
subj=adg$subj
# see the relationship between avd and intwt which are jointly modeled
library(lattice)
xyplot(y[,1]~log(y[,2]) | subj, ylab="Average Daily Gain",xlab="Initial Weight")
# below adg$subj is the block or barn
subj<-adg$subj
pred <- cbind(adg$pred.int,adg$pred.dummy1,adg$pred.dummy2,adg$pred.dummy3)
xcol<-1:4
zcol<-1
unst.psi.result <- mlmmm.em(y,subj,pred,xcol,zcol,maxits=200,eps=0.0001)
}}

\author{
Recai M. Yucel, Division of Biostatistics and Epidemiology, University of Massachusetts-Amherst
\email{yucel@schoolph.umass.edu}.}
\keyword{models}

% Converted by Sd2Rd version 37351.
