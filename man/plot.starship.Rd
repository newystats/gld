\name{plot.starship}
\alias{plot.starship}
\title{Plots to compare a fitted generalised lambda distribution to data}
\description{Plots to compare a Generalised Lambda Distribution fitted via 
the \code{\link{starship}} to data}
\usage{
\method{plot}{starship}(x, data, ask = FALSE, one.page = TRUE, 
    breaks = "Sturges", histogram.title = NULL,...)
}
\arguments{
  \item{x}{An object of class  \code{\link{starship}}.}
  \item{data}{Data to which the gld was fitted.  Leave this as NULL if the return.data argument was TRUE in the \code{\link{starship}} call that created x.}
  \item{ask}{Ask for user input before next plot --- passed to \code{\link{par}(ask)}.  Does not permanently change this setting. Ignored if \code{one.page} is TRUE}
  \item{one.page}{Put the two plots on one page using \code{\link{par}(mfrow=c(2,1))}. Does not permanently change this setting.}
  \item{breaks}{Control the number of histogram bins --- passed to \code{\link{hist}}.}
  \item{histogram.title}{Main title for histogram --- passed to \code{\link{hist}(main=)}.}
  \item{\dots}{arguments passed to \code{\link{plot}} AND \code{\link{hist}}}
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

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}
}
\seealso{
\code{\link{starship}},
}
\examples{
data <- rgl(100,0,1,.2,.2)
starship.result <- starship(data,optim.method="Nelder-Mead",initgrid=list(lcvect=(0:4)/10,
ldvect=(0:4)/10),return.data=TRUE)
plot(starship.result)
}
\keyword{distribution}


