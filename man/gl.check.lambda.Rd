\name{gl.check.lambda}
\alias{gl.check.lambda}
\title{Function to check the validity of parameters of the generalized lambda 
distribution}
\description{
  Checks the validity of parameters of the generalized lambda.  The tests 
  are simple for the FMKL, FM5 and GPD types, and  much more 
  complex for the RS parameterisation.
}
\usage{
gl.check.lambda(lambdas, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, param = "fkml", 
  lambda5 = NULL, vect = FALSE)
}
\arguments{
  \item{lambdas}{This can be either a single numeric value or a vector.  
  
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
  \item{lambda2}{\eqn{\lambda_2}{lambda 2} - scale parameter (\eqn{\beta} for \code{gpd})} 
  \item{lambda3}{\eqn{\lambda_3}{lambda 3} - first shape parameter (\eqn{\delta}, skewness parameter for \code{gpd})}
  \item{lambda4}{\eqn{\lambda_4}{lambda 4} - second shape parameter  (\eqn{\lambda}, kurtosis parameter for \code{gpd})}
  \item{lambda5}{\eqn{\lambda_5}{lambda 5} - a skewing parameter, in the
    fm5 parameterisation}
  \item{param}{choose parameterisation:
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}
  \code{fm5} uses the 5 parameter version of the FMKL parameterisation 
  (paper to appear)}
  \item{vect}{A logical, set this to TRUE if the parameters are given in the
  vector form (it turns off checking of the format of \code{lambdas} and the 
  other lambda arguments}
}
\details{
  See \code{\link{GeneralisedLambdaDistribution}} for details on the 
  generalised lambda distribution.  This function determines the validity of 
  parameters of the distribution.
  
  The FMKL parameterisation gives a valid 
  statistical distribution for any real values of \eqn{\lambda_1}{lambda 1},
  \eqn{\lambda_3}{lambda 3},\eqn{\lambda_4}{lambda 4} and any positive real 
  values of \eqn{\lambda_2}{lambda 2}.
  
  The FM5 parameterisation gives statistical distribution for any real 
  values of \eqn{\lambda_1}{lambda 1}, \eqn{\lambda_3}{lambda 3}, 
  \eqn{\lambda_4}{lambda 4}, any positive real values of 
  \eqn{\lambda_2}{lambda 2} and values of \eqn{\lambda_5}{lambda 5} that 
  satisfy \eqn{-1 \leq \lambda_5 \leq 1}{-1 <= lambda5 <= 1}.

  For the RS parameterisation, the combinations of parameters value that give 
  valid distributions are the following (the region numbers in the table 
  correspond to the labelling of the regions in \emph{Ramberg and 
  Schmeiser (1974)} and \emph{Karian, Dudewicz and McDonald (1996)}):

\tabular{cccccc}{
region \tab \eqn{\lambda_1}{lambda 1} \tab \eqn{\lambda_2}{lambda 2}	
\tab \eqn{\lambda_3}{lambda 3} \tab \eqn{\lambda_4}{lambda 4} \tab note \cr
1 \tab	all 	\tab \eqn{<0}	\tab \eqn{< -1}		\tab \eqn{> 1} \tab \cr
2 \tab	all 	\tab \eqn{<0}	\tab \eqn{> 1}		\tab \eqn{< -1} \tab \cr
3 \tab	all 	\tab \eqn{>0}	\tab \eqn{\ge 0}	\tab \eqn{\ge 0} \tab 
one of \eqn{\lambda_3}{lambda 3} and \eqn{\lambda_4}{lambda 4} must be non-zero \cr
4 \tab  all     \tab \eqn{<0} 	\tab \eqn{\le 0}	\tab \eqn{\le 0} \tab 
one of \eqn{\lambda_3}{lambda 3} and \eqn{\lambda_4}{lambda 4} must be non-zero \cr
5 \tab 	all	\tab \eqn{<0}	\tab \eqn{> -1} and \eqn{< 0} 	\tab \eqn{>1}
\tab equation 1 below must also be satisfied \cr
6 \tab 	all	\tab \eqn{<0}	\tab \eqn{>1} 	\tab \eqn{> -1} and \eqn{< 0}
\tab equation 2 below must also be satisfied \cr
}

Equation 1

\deqn{  \frac{(1-\lambda_3) ^{1-\lambda_3}(\lambda_4-1)^{\lambda_4-1}}
{(\lambda_4-\lambda_3)^{\lambda_4-\lambda_3}} <
- \frac{\lambda_3}{\lambda_4} }{% non TeX version
( (1-lambda3) ^ ( 1 - lambda3) * (lambda4 -1) ^ (lambda4 -1) ) /
( (lambda4 - lambda3) ^ (lambda4 - lambda3) ) <
- lambda3 / lambda 4 }

Equation 2

\deqn{  \frac{(1-\lambda_4) ^{1-\lambda_4}(\lambda_3-1)^{\lambda_3-1}}
{(\lambda_3-\lambda_4)^{\lambda_3-\lambda_4}} <
- \frac{\lambda_4}{\lambda_3} }{%
( (1-lambda4) ^ ( 1 - lambda4) * (lambda3 -1) ^ (lambda3 -1) ) /
( (lambda3 - lambda4) ^ (lambda3 - lambda4) ) <
- lambda4 / lambda 3
}

  The GPD type gives a valid distribution provided \eqn{\beta} is 
  positive and \eqn{0 \leq \delta \leq 1}{0 <= delta <= 1}.
}
\value{This logical function takes on a value of TRUE if the parameter values 
given produce a valid statistical distribution and FALSE if they don't}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988),
        \emph{A study of the generalized tukey lambda family}, Communications
       in Statistics - Theory and Methods \bold{17}, 3547--3567.

Karian, Z.E., Dudewicz, E.J., and McDonald, P. (1996), \emph{The extended 
    generalized lambda distribution system for fitting distributions to data:
    history, completion of theory, tables, applications, the ``Final Word'' 
    on Moment fits}, Communications in Statistics - Simulation and
Computation 
    \bold{25}, 611--642.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17},
78--82.

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{http://tolstoy.newcastle.edu.au/~rking/}}
\note{The complex nature of the rules in this function for the RS
parameterisation are the reason for the invention of the FMKL
parameterisation and its status as the default parameterisation in the other
generalized lambda functions.}

\seealso{The generalized lambda functions \code{\link{GeneralisedLambdaDistribution}}  }

\examples{
gl.check.lambda(c(0,1,.23,4.5),vect=TRUE) ## TRUE
gl.check.lambda(c(0,-1,.23,4.5),vect=TRUE) ## FALSE 
gl.check.lambda(c(0,1,0.5,-0.5),param="rs",vect=TRUE) ## FALSE
gl.check.lambda(c(0,2,1,3.4,1.2),param="fm5",vect=TRUE) ## FALSE
}
\keyword{distribution}
