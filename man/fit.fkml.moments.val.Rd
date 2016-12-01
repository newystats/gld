\name{fit.fkml.moments.val}
\alias{fit.fkml.moments.val}
\title{Method of moments estimation for the
FKML type of the generalised lambda distribution using given moment values} 
\description{Estimates parameters of the generalised lambda distribution
(FKML type) using the Method of Moments on the basis of moment values 
(mean, variance, skewness ratio and kurtosis ratio (note, not the \emph{excess kurtosis})). 
}
\usage{
fit.fkml.moments.val(moments=c(mean=0, variance=1, skewness=0,
  kurtosis=3), optim.method="Nelder-Mead", optim.control= list(), 
  starting.point = c(0,0))
}
\arguments{
  \item{moments}{A vector of length 4, consisting of the mean, variance
  and moment ratios for skewness and kurtosis (do not subtract 3 from the kurtosis ratio)}
  \item{optim.method}{Optimisation method for \code{\link{optim}} to use, 
defaults to Nelder-Mead}
  \item{optim.control}{argument \code{control}, passed to \code{\link{optim}}.}
  \item{starting.point}{a vector of length 2, giving the starting value
  for \eqn{\lambda_3}{lambda 3} and \eqn{\lambda_4}{lambda 4}.}
}
\details{Estimates parameters of the generalised lambda distribution
(FKML type) using Method of Moments on the basis of moment values (mean,
variance, skewness ratio and kurtosis ratio).  Note this is the fourth
central moment divided by the second central moment, without subtracting 3.
\code{fit.fkml.moments} (to come in version 2.4 of the \code{gld} package) will estimate
using the method of moments for a dataset, including calculating the
sample moments.
This function uses \code{\link{optim}} to find the parameters that
minimise the sum of squared differences between the skewness and kurtosis
sample ratios and their counterpart expressions for those ratios on the
basis of the parameters \eqn{\lambda_3}{lambda 3} and
\eqn{\lambda_4}{lambda 4}.
On the basis of these estimates (and the mean and variance), this function
then estimates 
\eqn{\hat\lambda_2}{lambda 2} and then \eqn{\hat\lambda_1}{lambda 1}.

Note that the first 4 moments don't uniquely identify members of the 
generalised \eqn{\lambda}{lambda} distribution.  Typically,
for a set of moments that correspond to a unimodal gld, there is
another set of parameters that give a distrbution with the same first 
4 moments.  This other distribution has a truncated appearance 
(that is, the distribution has finite support and the density is non-zero 
at the end points).  See the examples below.
}
\value{
A vector containing the parameters of the FKML type generalised lambda;
  \eqn{\lambda_1}{lambda 1} - location parameter
  \eqn{\lambda_2}{lambda 2} - scale parameter
  \eqn{\lambda_3}{lambda 3} - first shape parameter
  \eqn{\lambda_4}{lambda 4} - second shape parameter
(See \code{\link[gld]{gld}} for more details)
}
\references{
Au-Yeung, Susanna W. M. (2003) \emph{Finding Probability Distributions From 
Moments}, Masters thesis, Imperial College of Science, Technology and Medicine (University of London), Department of Computing
\url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.106.6130&rep=rep1&type=pdf}

Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
  \emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.
  
Lakhany, Asif and Mausser, Helmut (2000)
\emph{Estimating the parameters of the generalized lambda distribution},
Algo Research Quarterly, \bold{3(3)}:47--58

van Staden, Paul (2013) \emph{Modeling of generalized families of 
probability distributions inthe quantile statistical universe},
PhD thesis, University of Pretoria.
\url{http://repository.up.ac.za/handle/2263/40265}

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}

Sigbert Klinke

Paul van Staden}
\seealso{ 
\code{\link{gld.moments}} 
}
\examples{
# Approximation to the standard normal distribution
norm.approx <- fit.fkml.moments.val(c(0,1,0,3))
norm.approx
# Another distribution with the same moments
another <- fit.fkml.moments.val(c(0,1,0,3),start=c(2,2))
# Compared
plotgld(norm.approx,ylim=c(0,0.75),main="Approximation to the standard normal",
  sub="and another GLD with the same first 4 moments")
plotgld(another,add=TRUE,col=2)
}
\keyword{distribution}


