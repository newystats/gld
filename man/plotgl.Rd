\name{plotgl}
\alias{plotgl}
\alias{plotglc}
\alias{plotgld}
\title{Plots of density and distribution function for the generalised lambda
distribution}
\description{
Produces plots of density and distribution function for the generalised lambda
distribution.  Although you could use \code{plot(function(x)dgl(x))} to do
this, the fact that the density and quantiles of the generalised lambda are
defined in terms of the depth, \eqn{u}, means that a seperate function that
uses the depths to produce the values to plot is more efficient
}
\usage{
plotgld(lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fmkl", lambda5 = NULL, add = NULL, truncate = 0, 
  bnw = FALSE, col.or.type = 1, granularity = 10000, xlab = "x", 
  ylab = NULL, quant.probs = seq(0,1,.25), new.plot = NULL,  ...)
plotglc(lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fmkl", lambda5 = NULL, granularity = 10000, xlab = "x", 
  ylab = "cumulative probability", add = FALSE, ...)
}
\arguments{
  \item{lambda1}{This can be either a single numeric value or a vector.  
  
  If it is a vector, it must be of length 4 for parameterisations 
  \code{fmkl} or \code{rs} and of length 5 for parameterisation \code{fm5}.  
  If it is a vector, it gives all the parameters of the generalised lambda 
  distribution (see below for details) and the other \code{lambda} arguments
  must be left as NULL.
        
 If it is a a single value, it is \eqn{\lambda_1}{lambda 1}, the location 
 parameter of the distribution and the other parameters are given by the
 following arguments
 
 \emph{Note that the numbering of the \eqn{\lambda}{lambda} parameters for
 the \emph{fmkl} parameterisation is different to that used by Freimer, 
 Mudholkar, Kollia and Lin.}
	}
  \item{lambda2}{\eqn{\lambda_2}{lambda 2} - scale parameter}
  \item{lambda3}{\eqn{\lambda_3}{lambda 3} - first shape parameter}
  \item{lambda4}{\eqn{\lambda_4}{lambda 4} - second shape parameter}
  \item{lambda5}{\eqn{\lambda_5}{lambda 5} - a skewing parameter, in the
    fm5 parameterisation}
  \item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}
  \code{fm5} uses the 5 parameter version of the FMKL parameterisation 
  (paper to appear)}
  \item{add}{a logical value describing whether this should add to an existing plot (using
\code{lines}) or produce a new plot (using \code{plot}). Defaults to FALSE (new plot) if
both \code{add} and \code{new.plot} are NULL.}
  \item{truncate}{for \code{plotgld}, a minimum density value at which the
plot should be truncated.}
  \item{bnw}{a logical value, true for a black and white plot, with different
densities identified using line type (\code{lty}), false for a colour plot, 
with different
densities identified using line colour (\code{col})}
  \item{col.or.type}{Colour or type of line to use}
  \item{granularity}{Number of points to calculate quantiles and density at
--- see \emph{details}}
  \item{xlab}{X axis label}
  \item{ylab}{Y axis label}
  \item{quant.probs}{Quantiles of distribution to return (see \emph{value}
below).  Set to NULL to suppress this return entirely.}
  \item{new.plot}{a logical value describing whether this should produce a
new plot (using \code{plot}), or add to an existing plot (using
\code{lines}).  Ignored if \code{add} is set.}
  \item{...}{arguments that get passed to \code{plot} if this is a new plot}
}
\details{
The generalised lambda distribution is defined in terms of its quantile
function.  The density of the distribution is available explicitly as a
function of depths, \eqn{u}, but not explicitly available as a function of 
\eqn{x}.  This function calculates quantiles and depths as a function of
depths to produce a density plot \code{plotgld} or cumulative probability plot
\code{plotglc}.

The plot can be truncated, either by restricting the values using \code{xlim}
--- see \code{par} for details, or by the \code{truncate} argument, which
specifies a minimum density.  This is recommended for graphs of densities
where the tail is very long.
}
\value{
  A number of quantiles from the distribution, the default being the minimum,
maximum and quartiles.
}
\references{ 
 Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.

Karian, Z.E. & Dudewicz, E.J. (2000), \emph{Fitting Statistical Distributions to
Data: The generalised Lambda Distribution and the Generalised Bootstrap
Methods}, CRC Press.

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}

}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}}

\seealso{\code{\link{GeneralisedLambdaDistribution}}}
\examples{
plotgld(0,1.4640474,.1349,.1349,main="Approximation to Standard Normal", 
sub="But you can see this isn't on infinite support")

plotgld(1.42857143,1,.7,.3,main="The whale")
plotglc(1.42857143,1,.7,.3)
plotgld(0,-1,5,-0.3,param="rs")
plotgld(0,-1,5,-0.3,param="rs",xlim=c(1,2))
# A bizarre shape from the RS paramterisation
plotgld(0,1,5,-0.3,param="fmkl")
plotgld(10/3,1,.3,-1,truncate=1e-3)

plotgld(0,1,.0742,.0742,col.or.type=2,param="rs",
main="All distributions have the same moments",
sub="The full Range of all distributions is shown")
plotgld(0,1,6.026,6.026,col.or.type=3,new.plot=FALSE,param="rs")
plotgld(0,1,35.498,2.297,col.or.type=4,new.plot=FALSE,param="rs")
legend(0.25,3.5,lty=1,col=c(2,3,4),legend=c("(0,1,.0742,.0742)",
"(0,1,6.026,6.026)","(0,1,35.498,2.297)"),cex=0.9)
# An illustration of problems with moments as a method of characterising shape

}
\keyword{distribution}
\keyword{hplot}
\keyword{aplot}

