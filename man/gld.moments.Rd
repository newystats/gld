\name{gld.moments}
\alias{gld.moments} 
\title{Calculate moments of the FKML type of the generalised lambda
distribution for given parameter values} 
\description{Calculates the mean, variance, skewness ratio and kurtosis ratio of the generalised \eqn{\lambda}{lambda} distribution for
given parameter values.
}
\usage{
  gld.moments(par,type="fkml",ratios=TRUE)
}
\arguments{
  \item{par}{A vector of length 4, giving the parameters of the 
  generalised lambda distribution, consisting of; 
  {\eqn{\lambda_1}{lambda 1} location parameter}
  {\eqn{\lambda_2}{lambda 2} - scale parameter}
  {\eqn{\lambda_3}{lambda 3} - first shape parameter}
  {\eqn{\lambda_4}{lambda 4} - second shape parameter}
    }
  \item{type}{choose the type of generalised lambda distribution.  Currently  \code{gld.moments} only supports
  \code{fkml} which uses \emph{Freimer, Kollia, Mudholkar, and Lin (1988)} (default).}
  \item{ratios}{Logical.  TRUE to give moment ratios for skewness and kurtosis, FALSE to give the third and fourth central moments instead.}
}
\details{
The FKML type of the generalised \eqn{\lambda}{lambda} distribution was
introduced by Freimer et al (1988) who gave expressions for the moments.
In the limit, as the shape parameters (\eqn{\lambda_3}{lambda 3} and 
\eqn{\lambda_4}{lambda 4}) go to zero, the distribution is defined using
limit results.  The moments in these limiting cases were given by 
van Staden (2013).  This function calculates the first 4 moments.

See pages 96--97 of van Staden (2013) for the full expressions for these
moments.
}
\value{
A vector containing the first four moments of the FKML type generalized
lambda.  If \code{ratio} is true, the vector contains the mean,
variance, skewness ratio and kurtosis ratio.  If \code{ratio} is false,
the vector contains the mean, variance, third central moment and fourth
central moment.
}
\references{
Au-Yeung, Susanna W. M. (2003) \emph{Finding Probability Distributions From 
Moments}, Masters thesis, Imperial College of Science, Technology and Medicine (University of London), Department of Computing
\url{https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.106.6130&rep=rep1&type=pdf}

Freimer, M., Kollia, G., Mudholkar, G. S.,  & Lin, C. T. (1988), 
  \emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.
  
Lakhany, Asif and Mausser, Helmut (2000)
\emph{Estimating the parameters of the generalized lambda distribution},
Algo Research Quarterly, \bold{3(3)}:47--58

van Staden, Paul J. (2013) \emph{Modeling of generalized families of 
probability distributions inthe quantile statistical universe},
PhD thesis, University of Pretoria.
\url{https://repository.up.ac.za/handle/2263/40265}

\url{https://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{https://tolstoy.newcastle.edu.au/~rking/}

Sigbert Klinke

Paul van Staden}
\seealso{\code{\link{fit.fkml.moments.val}}
}
\examples{
gld.moments(c(0,1.463551,0.1349124,0.1349124))
gld.moments(c(0,1.813799,0,0))
gld.moments(c(0,1,0,3))
}
\keyword{distribution}


