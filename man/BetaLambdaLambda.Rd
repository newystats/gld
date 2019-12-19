\name{BetaLambdaLambda}
\alias{BetaLambdaLambda}
\title{Calculates Beta function for two identical parameters, allowing non-integer negative values}
\description{By defining the Beta Function in terms of the Gamma
Function,
\deqn{B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}}{
Beta(a,b)=(Gamma(a)*Gamma(b))/Gamma(a+b)}
the function can be defined for non-integer negative values of 
a and b.  The special case of this where \eqn{a=b} is needed 
to calculate the standard errors of the L Moment estimates of the
gpd type of the generalised lambda distribution, so this function
carries out that calculation.
}
\usage{
BetaLambdaLambda(lambda)
}
\arguments{
  \item{lambda}{A vector, each element of which is used for both 
  arguments of the Beta function.}
}
\details{\code{NaN} is returned for any negative integer elements of \code{lambda}.
}
\value{
A vector the same length as \code{lambda}, containing 
  Beta(lambda,lambda)
}
\references{
\url{https://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{https://tolstoy.newcastle.edu.au/~rking/}

Paul van Staden
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{beta}}
\code{\link{gamma}}
\code{\link{fit.gpd}} 
\code{\link{GeneralisedLambdaDistribution}}
}
\examples{
BetaLambdaLambda(-0.3)
}
\keyword{math}
