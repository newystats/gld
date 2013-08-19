\name{GeneralisedLambdaDistribution}
\alias{dgl}
\alias{qdgl}
\alias{dqgl}
\alias{pgl}
\alias{qgl}
\alias{rgl}
\alias{GeneralizedLambdaDistribution}
\alias{GeneralisedLambdaDistribution}
\alias{gld}
\title{The Generalised Lambda Distribution}
\description{Density, quantile density, distribution function, quantile 
function and random generation for the generalised lambda distribution 
(also known as the asymmetric lambda, or Tukey lambda).  Provides for three
different parameterisations, the \code{fmkl} (recommended), the \code{rs} and 
a five parameter version of the FMKL, the \code{fm5}.
}
\usage{
dgl(x, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fkml", lambda5 = NULL, inverse.eps = .Machine$double.eps,
  max.iterations = 500)
dqgl(p, lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fkml", lambda5 = NULL)
qdgl(p, lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fkml", lambda5 = NULL)
pgl(q, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
  param = "fkml", lambda5 = NULL, inverse.eps = .Machine$double.eps,
  max.iterations = 500)
qgl(p, lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL,
  param = "fkml", lambda5 = NULL)
rgl(n, lambda1=0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL,
  param = "fkml", lambda5 = NULL)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
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
  \item{inverse.eps}{Accuracy of calculation for the numerical determination of 
	\eqn{F(x)}, defaults to \code{.Machine$double.eps}}
  \item{max.iterations}{Maximum number of iterations in the numerical 
	determination of \eqn{F(x)}, defaults to 500}
}
\details{
The generalised lambda distribution, also known as the asymmetric lambda, 
or Tukey lambda distribution, is a distribution with a wide range of shapes.  
The distribution is defined by its quantile function, the inverse of the
distribution function.  The \code{gld} package implements three parameterisations of the distribution.  
The default parameterisation (the FMKL) is that due to \emph{Freimer 
Mudholkar, Kollia and Lin (1988)} (see references below), with a quantile 
function:
\deqn{F^{-1}(u)= \lambda_1 + { { \frac{u^{\lambda_3}-1}{\lambda_3} - 	%
\frac{(1-u)^{\lambda_4}-1}{\lambda_4} } \over \lambda_2 } 		%
}{F inverse (u) = lambda1 + ( (u^lambda3 -1)/lambda3 - ((1-u)^lambda4-1) 
/lambda4 ) / lambda 2}
for \eqn{\lambda_2 > 0}{lambda2 >0}.

A second parameterisation, the RS, chosen by setting \code{param="rs"} is 
that due to \emph{Ramberg and Schmeiser (1974)}, with the quantile function:
\deqn{F^{-1}(u)= \lambda_1 + \frac{u^{\lambda_3} - (1-u)^{\lambda_4}} 	%
{\lambda_2} }{F inverse (u) = lambda1 + ( u^lambda3 - (1-u)^lambda4 ) / % 
lambda 2}

This parameterisation has a complex series of rules determining which values 
of the parameters produce valid statistical distributions.  See 
\link{gl.check.lambda} for details.

A third parameterisation, the FM5, chosen by setting \code{param="fm5"} 
adds an additional skewing parameter to the FMKL parameterisation.  
This uses the same approach as that used by 
Gilchrist (2000)
for the RS parameterisation.  The quantile function is
\deqn{F^{-1}(u)= \lambda_1 + { { \frac{(1-\lambda_5)(u^{\lambda_3}-1)}%
{\lambda_3} - \frac{(1+\lambda_5)((1-u)^{\lambda_4}-1)}{\lambda_4} }  %
\over \lambda_2 }}{F inverse (u) = lambda1 + ( (1-lambda5)(u^lambda3-1)/lambda3 
- (1+lambda5)((1-u)^lambda4-1)/lambda4 ) / lambda 2}
for \eqn{\lambda_2 > 0}{lambda2 >0} 
and \eqn{-1 \leq \lambda_5 \leq 1}{-1 <= lambda5 <= 1}.

The distribution is defined by its quantile function and its distribution and
density functions do not exist in closed form.  Accordingly, the results 
from \code{pgl} and \code{dgl} are the result of numerical solutions to the
quantile function, using the Newton-Raphson method.  Since the quantile density
function, \eqn{f(F^{-1}(u))}, does exist, an additional function, \code{qdgl},
computes this.

The functions \code{qdgl.fmkl}, \code{qdgl.rs}, \code{qdgl.fm5},
\code{qgl.fmkl}, \code{qgl.rs} and \code{qgl.fm5} from versions 1.5 and 
earlier of the \code{gld} package have been renamed (and hidden) to \code{.qdgl.fmkl}, 
\code{.qdgl.rs}, \code{..qdgl.fm5}, \code{.qgl.fmkl}, \code{.qgl.rs} 
and \code{.qgl.fm5} respectively.  See the code for more details

}
\value{
  \code{dgl} gives the density (based on the quantile density and a 
  numerical solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}), 

  \code{qdgl} gives the quantile density,
  
  \code{pgl} gives the distribution function (based on a numerical
	solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}),

  \code{qgl} gives the quantile function, and

  \code{rgl} generates random deviates.
}
\references{
Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Gilchrist, Warren G. (2000), \emph{Statistical Modelling with Quantile Functions}, Chapman \& Hall

Karian, Z.A., Dudewicz, E.J., and McDonald, P. (1996), \emph{The extended 
    generalized lambda distribution system for fitting distributions to data: 
    history, completion of theory, tables, applications, the ``Final Word'' 
    on Moment fi
ts}, Communications in Statistics - Simulation and Computation 
    \bold{25}, 611--642.
    
Karian,  Zaven A. and Dudewicz, Edward J. (2000), \emph{Fitting statistical distributions: the Generalized Lambda Distribution and Generalized Bootstrap methods}, Chapman & Hall

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{http://tolstoy.newcastle.edu.au/~rking/}}

\examples{
qgl(seq(0,1,0.02),0,1,0.123,-4.3)
pgl(seq(-2,2,0.2),0,1,-.1,-.2,param="fmkl",inverse.eps=.Machine$double.eps)
# calculate the probabilities more accurately than normal
}
\keyword{distribution}
