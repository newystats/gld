\name{starship.adaptivegrid}
\alias{starship.adaptivegrid}
\title{Carry out the ``starship'' estimation method for the generalised
lambda distribution using a grid-based search}
\description{
Calculates estimates for the generalised lambda
distribution on the basis of data, using the starship method.
The starship method is built on the fact that the 
generalised lambda distribution (\code{\link{gld}})
is a transformation of the uniform distribution.  This method finds the 
parameters that transform the data closest to the uniform distribution.
This function uses a grid-based search.
}
\usage{
starship.adaptivegrid(data, initgrid,inverse.eps = 1e-08, param="FMKL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data to be fitted, as a vector}
  \item{initgrid}{A list with elements, 
  \code{lcvect}, a vector of values for \eqn{\lambda_3}{lambda 3},
  \code{ldvect}, a vector of values for \eqn{\lambda_4}{lambda 4} and
  \code{levect}, a vector of values for \eqn{\lambda_5}{lambda 5} 
  (\code{levect} is only required if \code{param} is \code{fm5}).
  The parameter values given in \code{initgrid} are not checked with \code{\link{gl.check.lambda}}.
  }
  \item{inverse.eps}{Accuracy of calculation for the numerical determination of 
	\eqn{F(x)}, defaults to \eqn{10^{-8}}{1e-8}}
  \item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}
  \code{fm5} uses the 5 parameter version of the FMKL parameterisation 
  (paper to appear)}
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

This function carries out a grid-based search.  This was the original method
of King \& MacGillivray, 1999, but you are advised to instead use 
\code{\link{starship}} which uses a grid-based search together with an
optimisation based search.

See \code{\link{GeneralisedLambdaDistribution}} for details on
parameterisations.
}
\value{
  \item{response}{The minimum ``response value'' --- the result of the
internal goodness-of-fit measure.  This is the return value of starship.obj.
See King \& MacGillivray, 1999 for more details}
  \item{lambda}{A vector of length 4 giving the values of
\eqn{\lambda_1}{lambda 1} to \eqn{\lambda_4}{lambda 4} that produce this
minimum response, i.e. the estimates}
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

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}

Darren Wraith}
\seealso{ 
\code{\link{starship}}, 
\code{\link{starship.obj}} 
}
\examples{
data <- rgl(100,0,1,.2,.2)
starship.adaptivegrid(data,list(lcvect=(0:4)/10,ldvect=(0:4)/10))
}
\keyword{distribution}
