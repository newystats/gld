\name{starship.obj}
\alias{starship.obj}
\title{Objective function that is minimised in starship estimation method}
\description{
The starship is a method for fitting the generalised lambda distribution.
See \code{\link{starship}} for more details.

This function is the objective funciton minimised in the methods.  It is a
goodness of fit measure carried out on the depths of the data.
}
\usage{
starship.obj(par, data, inverse.eps, param = "fmkl")
}
\arguments{
  \item{par}{parameters of the generalised lambda distribution, a vector of
length 4, giving \eqn{\lambda_1}{lambda 1} to \eqn{\lambda_4}{lambda 4}.
See \code{\link{GeneralisedLambdaDistribution}} for details on the
definitions of these parameters}
  \item{data}{Data --- a vector}
\item{inverse.eps}{Accuracy of calculation for the numerical determination of 
	\eqn{F(x)}, defaults to \eqn{10^{-8}}{1e-8}}
\item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}}
}
\details{
The starship method is described in King \& MacGillivray, 1999 (see
references). It is built on the fact that the 
generalised lambda distribution (\code{\link{gld}})
is a transformation of the uniform distribution.  Thus the inverse of this
transformation is the distribution function for the gld.  The starship method
applies different values of the parameters of the distribution to the
distribution function, calculates the depths \emph{q} corresponding to the data
and chooses the parameters that make the depths closest to a uniform
distribution.

The closeness to the uniform is assessed by calculating the Anderson-Darling
goodness-of-fit test on the transformed data against the uniform, for a
sample of size \code{length(data)}.

This function returns that objective function.  It is provided as a seperate
function to allow users to carry out minimisations using \code{\link{optim}}
or other methods.  The recommended method is to use the \code{\link{starship}}
function.
}
\value{
The Anderson-Darling goodness of fit measure, computed on the transformed
data, compared to a uniform distribution.  \emph{Note that this is NOT the
goodness-of-fit measure of the generalised lambda distribution
with the given parameter values to the data.}
}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.
 
King, R.A.R. & MacGillivray, H. L. (1999), \emph{A starship method for
fitting the generalised \eqn{\lambda}{lambda} distributions},
Australian and New Zealand Journal of
Statistics \bold{41}, 353--374

Owen, D. B. (1988), \emph{The starship}, Communications in Statistics - 
Computation and Simulation \bold{17}, 315--323.

\url{https://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{https://tolstoy.newcastle.edu.au/~rking/}
Darren Wraith}
\seealso{ 
\code{\link{starship}},
\code{\link{starship.adaptivegrid}}
}
\examples{
data <- rgl(100,0,1,.2,.2)
starship.obj(c(0,1,.2,.2),data,inverse.eps=1e-10,"fmkl")
}
\keyword{distribution}
