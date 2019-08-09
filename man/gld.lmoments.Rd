\name{gld.lmoments}
\alias{gld.lmoments} 
\title{Calculate L-Moments of the GPD type generalised lambda
distribution for given parameter values} 
\description{Calculates the first four L-Moments of the GPD type 
generalised \eqn{\lambda}{lambda} distribution for given parameter 
values.}
\usage{
  gld.lmoments(pars,order=1:4,ratios=TRUE,type="GPD",param=NULL)
}
\arguments{
  \item{pars}{A vector of length 4, giving the parameters of the 
  GPD type generalised lambda distribution, consisting of; 
  \itemize{
    \item{\eqn{\alpha}{alpha}}{ location parameter}
    \item{\eqn{\beta>0}{beta > 0}}{ scale parameter}
    \item{\eqn{0 \leq \delta \leq 1}{0 <= delta <= 1}}{ skewness parameter}
    \item{\eqn{\lambda}{lambda}}{ kurtosis parameter}
      }
    }
  \item{order}{Integers to select the orders of L-moments to calculate.
    Currently this function only calculates for orders 1 to 4.}
  \item{type}{choose the type of generalised lambda distribution.  Currently  \code{gld.lmoments} only supports
  \code{GPD} which uses \emph{van Staden and Loots (2009)} (default).}
  \item{ratios}{Logical.  TRUE gives L-moment ratios for skewness and 
    kurtosis (\eqn{\tau_3}{tau 3} and \eqn{\tau_4}{tau 4}) (and all higher
    orders), 
    FALSE gives the requested L-moments instead.}
  \item{param}{alias for the type argument.  The \code{type} argument is preferred.}
}
\details{
The GPD type generalised \eqn{\lambda}{lambda} distribution was
introduced by van Staden and Loots (2009).  It has explicit parameters
for skewness and kurtosis, and closed form estimates for L-moment 
estimates of the parameters.

In the limit, as the kurtosis parameter, \eqn{\lambda}{lambda}, goes 
to zero, the distribution approaches the skew logistic distribution of
van Staden and King (2013).  See the \code{sld} package for this 
distribution.
}
\value{
A vector containing the selected L-moments of the GPD type generalised
lambda.  If \code{ratio} is true, the vector contains L-Moment ratios for orders 3 and over, otherwise all values are L-Moments.
}
\references{
van Staden, P.J. and King, Robert A.R. (2015) \emph{The quantile-based skew logistic distribution}, Statistics and Probability Letters \bold{96} 109--116.
\url{http://dx.doi.org/10.1016/j.spl.2014.09.001}

van Staden, Paul J. 2013 \emph{Modeling of generalized families of probability distribution in the quantile statistical universe}. 
PhD thesis, University of Pretoria.
\url{http://hdl.handle.net/2263/40265}

Van Staden, Paul J., & M.T. Loots. (2009), \emph{Method of L-moment Estimation for the Generalized Lambda Distribution}. In Proceedings of the Third Annual ASEARC Conference. Callaghan, NSW 2308 Australia: School of Mathematical and Physical Sciences, University of Newcastle.

Quantile based Skew logistic distribution

\url{http://tolstoy.newcastle.edu.au/rking/SLD/SLD.html}

Generalised Lambda Distribution

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}

Paul van Staden}
\seealso{
\code{sld} package
}
\examples{
gld.lmoments(c(0,1,0.5,0.23))
gld.lmoments(c(0,1,0,0.23))
gld.lmoments(c(0,1,0.5,0.7))
}
\keyword{distribution}


