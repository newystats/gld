\name{print.starship}
\alias{summary.starship}
\alias{print.starship}
\title{Print (or summarise) the results of a starship estimation}
\description{Print (or summarise) the results of a \code{\link{starship}} estimation
of the parameters of the Generalised Lambda Distribution}
\usage{
\method{summary}{starship}(object, ...)

\method{print}{starship}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{An object of class  \code{\link{starship}}.}
  \item{object}{An object of class  \code{\link{starship}}.}
  \item{digits}{minimal number of \emph{significant} digits, see
    \code{\link{print.default}}.}
  \item{\dots}{arguments passed to \code{\link{print}}}
}
\details{
\code{summary} Gives the details of the \code{\link{starship.adaptivegrid}} and \code{optim} steps.
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
\code{\link{starship.adaptivegrid}},
\code{\link{starship.obj}} 
}
\examples{
data <- rgl(100,0,1,.2,.2)
starship.result <- starship(data,optim.method="Nelder-Mead",initgrid=list(lcvect=(0:4)/10,
ldvect=(0:4)/10))
print(starship.result)
summary(starship.result,estimation.details=TRUE)
}
\keyword{distribution}

