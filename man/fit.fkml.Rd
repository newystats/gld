\name{fit.fkml}
\alias{fit.fkml}
\title{Estimate parameters of the FKML parameterisation of the generalised lambda distribution
}
\description{
Estimates parameters of the FKML parameterisation of the Generalised \eqn{\lambda}{lambda} Distribution.  Five estimation methods are available;
Numerical Maximum Likelihood, 
Maximum Product of Spacings, Titterington's Method, 
the Starship (also available in the \code{\link{starship}} function, which uses the same underlying code as this for the fkml parameterisation),
and Trimmed L-Moments.
}
\usage{
fit.fkml(x, method = "ML", t1 = 0, t2 = 0, 
  l3.grid = c(-0.9, -0.5, -0.1, 0, 0.1, 0.2, 0.4, 0.8, 1, 1.5), 
  l4.grid = l3.grid, record.cpu.time = TRUE, optim.method = "Nelder-Mead",
  inverse.eps = .Machine$double.eps, optim.control=list(maxit=10000), 
  optim.penalty=1e20, return.data=FALSE)
}
\arguments{
  \item{x}{Data to be fitted, as a vector}
  \item{method}{A character string, to select the estimation method.  One of: \code{ML} for Numerical Maximum Likelihood, \code{MPS} or \code{MSP}
  for Maximum Product of Spacings, \code{TM} for Titterington's Method, \code{SM} for  Starship Method, \code{TL} for  Method of TL-moments, or \code{DLA} for the method of distributional least absolutes.}
  \item{t1}{Number of observations to be trimmed from the left in the conceptual sample, \eqn{t_1}{t1} (A non-negative integer, only used by TL-moment estimation, see details section)}
  \item{t2}{Number of observations to be trimmed from the right in the conceptual sample, \eqn{t_2}{t2} (A non-negative integer, only used by TL-moment estimation, see details section). These two arguments are 
  restricted by \eqn{t_1+t_2 < n}{t1+t2 < n}, where n is the sample size}
  \item{l3.grid}{A vector of values to form the grid of values of \eqn{\lambda_3}{lambda3} used to find a starting point for the optimisation.}
  \item{l4.grid}{A vector of values to form the grid of values of \eqn{\lambda_4}{lambda4} used to find a starting point for the optimisation.}
  \item{record.cpu.time}{Boolean --- should the CPU time used in fitting be recorded in the fitted model object?}
  \item{optim.method}{Optimisation method, use any of the options available under \code{method} of \code{\link[stats]{optim}}.}
\item{inverse.eps}{Accuracy of calculation for the numerical determination of \eqn{F(x)}, defaults to \code{.Machine$double.eps}.}
\item{optim.control}{List of options for the optimisation step.  See
\code{\link[stats]{optim}} for details.
}
\item{optim.penalty}{The penalty to be added to the objective function
if parameter values are proposed outside the allowed region}
\item{return.data}{Logical: Should the function return the data (from the argument \code{data})?}
}
\details{
Maximum Likelihood Estimation of the generalised lambda distribution (\code{\link{gld}}) proceeds by calculating the density of the data for candidate values of the parameters.  Because the
gld is defined by its quantile function, the method first numerically 
obtains F(x) by inverting Q(u), then obtains the density for that observation.

Maximum Product of Spacings estimation (sometimes referred to as Maximum 
Spacing Estimation, or Maximum Spacings Product)  finds the 
parameter values that
maximise the product of the spacings (the difference between successive 
depths, \eqn{F_\theta(x_{(i+1)})-F_\theta(x_{(i)})}{F(x_(i+1);theta)-F(x_(i);theta)}, where \eqn{F_\theta(x)}{F(x;theta)} is the distribution
function for the candidate values of the parameters).  See % King & Dean 2013,
Dean (2013) and Cheng & Amin (1981) for details.

Titterington (1985) remarked that MPS effectively added 
an ``extra observation''; there are N data points in the original 
sample, but N + 1 spacings in the expression maximised in MPS.
Instead of using spacings between transformed data points, so method \code{TM}
uses  spacings between transformed, adjacently-averaged, data points.
The spacings are given by \eqn{D_i=F_\theta(z_{(i)})-F_\theta(z_{(i-1)})}{Di=F(z(i);theta)- F(z(i-1);theta)},
where \eqn{\alpha_1 = z_0 < z_1 < \dots < z_n = \alpha_2}{alpha1 = z0 < z1 < ... 
<zn = alpha2}
and \eqn{z_i = (x_{(i)} +x_{(i+1)})/2}{z_i = (x_i + x_(i+1) )/2} 
for i = 1,2, ... n-1
(where \eqn{\alpha_1}{alpha1} and \eqn{\alpha_2}{alpha2} 
are the lower and upper bounds on the support of the distribution). 
This reduces the number of spacings to \eqn{n}{n} and achieves concordance 
with the original sample size. See Titterington (1985) and 
Dean (2013) for details.

The starship is built on the fact that the \eqn{g\lambda d}{gld}
is a transformation of the uniform distribution.  Thus the inverse of 
this transformation is the distribution function for the gld.  The 
starship method applies different values of the parameters of the distribution to
the distribution function, calculates the depths \emph{q} corresponding to 
the data and chooses the parameters that make these calculated depths 
closest (as measured by the Anderson-Darling statistic) to a uniform distribution.
See King & MacGillivray (1999) for details.

TL-Moment estimation chooses the values of the parameters that minimise the 
difference between the sample Trimmed L-Moments of the data and the Trimmed 
L-Moments of the fitted distribution.  
TL-Moments are based on inflating the conceptual sample size used in the definition of L-Moments.  The \code{t1} and \code{t2} arguments to the
function define the extent of trimming of the conceptual sample.  Thus,
the default values of  \code{t1=0} and \code{t2=0} reduce the TL-Moment
method to L-Moment estimation.  \code{t1} and \code{t2} give the number of 
observations to be trimmed (from the left and right respectively) from 
the conceptual sample of size \eqn{n+t_1+t_2}{n+t1+t2}.  These two
arguments should be non-negative integers, and \eqn{t_1+t_2 < n}{t1+t2 < n}, where n is the sample size.
See Elamir and Seheult (2003) for more on 
TL-Moments in general, Asquith, (2007) for TL-Moments of the RS parameterisation of the gld and Dean (2013) for more details on TL-Moment estimation of the gld.

The method of distributional least absolutes (DLA) minimises the sum of absolute deviations between the order statistics and their medians (based on the candidate parameters). See Dean (2013) for more information.
}
\value{
\code{fit.fkml} returns an object of \code{\link[base]{class}} 
\code{"starship"} (regardless of the estimation method used).

\code{print} prints the estimated values of the parameters, while
\code{\link[gld]{summary.starship}} prints these by default, but can also provide
details of the estimation process (from the components \code{grid.results},
\code{data} and \code{optim} detailed below).

The value of \code{fit.fkml} is a list containing the
following components:

  \item{lambda}{A vector of length 4, giving 
the estimated parameters, in order,
  \eqn{\lambda_1}{lambda 1} - location parameter
  \eqn{\lambda_2}{lambda 2} - scale parameter
  \eqn{\lambda_3}{lambda 3} - first shape parameter
  \eqn{\lambda_4}{lambda 4} - second shape parameter
  }
  \item{grid.results}{output from the grid search}
  \item{optim}{output from the optim search,
\code{\link{optim}} for details}
  \item{cpu}{A vector showing the computing time used, returned if 
    \code{record.cpu.time} is TRUE}     
  \item{data}{The data, if \code{return.data} is TRUE}
}
\references{
Asquith, W. H. (2007), \emph{L-Moments and TL-Moments of the Generalized Lambda
Distribution}, Computational Statistics & Data Analysis, \bold{51}, 4484--4496.

Cheng, R.C.H. & Amin, N.A.K. (1981), \emph{Maximum Likelihood Estimation of Parameters in the Inverse Gaussian Distribution, with Unknown Origin}, 
Technometrics, \bold{23(3)}, 257--263.
\url{http://www.jstor.org/stable/1267789}

Dean, B. (2013) \emph{Improved Estimation and Regression Techniques with the Generalised Lambda Distribution}, PhD Thesis, University of Newcastle
\url{http://nova.newcastle.edu.au/vital/access/manager/Repository/uon:13503}

Elamir, E. A. H., and Seheult, A. H. (2003), \emph{Trimmed L-Moments}, Computational Statistics & Data Analysis, \bold{43}, 299--314.

King, R.A.R. & MacGillivray, H. L. (1999), \emph{A starship method for
fitting the generalised \eqn{\lambda}{lambda} distributions},
Australian and New Zealand Journal of
Statistics \bold{41}, 353--374.

Titterington, D. M. (1985), \emph{Comment on `Estimating Parameters in Continuous Univariate Distributions'}, Journal of the Royal Statistical Society, Series B, \bold{47}, 115--116.
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}

Ben Dean, University of Newcastle
\email{benjamin.dean@uon.edu.au}
}
\seealso{
\code{\link{starship}} 
\code{\link{GeneralisedLambdaDistribution}}
}
\examples{
example.data <- rgl(200,c(3,1,.4,-0.1),param="fkml")
example.fit <- fit.fkml(example.data,"MSP",return.data=TRUE)
print(example.fit)
summary(example.fit)
plot(example.fit,one.page=FALSE)
}
\keyword{distribution}



