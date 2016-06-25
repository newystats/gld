\name{fit.gpd}
\alias{fit.gpd}
\alias{fit.gpd.lmom}
\alias{fit.gpd.lmom.given} 
\title{Estimate parameters of the GPD type generalised lambda distribution}
\description{Estimates parameters of the GPD type generalised \eqn{\lambda}{lambda} Distribution.  Currently, only estimation via method of L moments is implemented.

The Method of L-Moments estimates for the GPD type are the only estimates for any
generalised lambda distribution type with closed form expressions.
}
\usage{
fit.gpd(x, method = "LM", na.rm = TRUE, record.cpu.time = TRUE, return.data = FALSE)
fit.gpd.lmom(data, na.rm = TRUE)
fit.gpd.lmom.given(lmoms, n = NULL)
}
\arguments{
  \item{x}{Data to be fitted, as a vector}
  \item{method}{A character string, to select the estimation method.  Only Method of
  L-Moments "LM" is implemented.}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{record.cpu.time}{Boolean: should the CPU time used in fitting be recorded in the fitted model object?}
  \item{return.data}{Boolean: Should the function return the data (from the argument \code{x})?}
  \item{data}{Data to be fitted, as a vector}
  \item{lmoms}{A numeric vector containing two L-moments and two L-moment ratios, in the order \eqn{l_1,l_2,t_3,t_4}{l_1, l_2, t_3, t_4}.}
  \item{n}{the sample size, defaults to NULL}
}
\details{The method of L-Moments equates sample L-Moments with expressions 
for the L-Moments of the GPD type GLD.  Closed form expressions exist to give these estimates.

If the L-Moment ratio, \eqn{\tau_4}{tau4} is less than the minimum 
value that \eqn{\tau_4}{tau4} can obtain for the GPD generalised lambda
distribution;

\deqn{\tau_4^{(min)}=\frac{12-5\sqrt{6}}{12+5\sqrt{6}} \approx -0.0102051,}{(12-5*sqrt(6))/(12+5*sqrt(6)) = approx -0.0102051,}

there is no possible estimate, and this function returns

For many values there are two possible estimates, in two regions of the
GPD GLD parameter space, denoted region A and region B in van Staden (2013).

See page 154

This function returns both estimates.  In most cases, the 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, method = "LM", na.rm = TRUE, record.cpu.time = TRUE, 
    return.data = FALSE) 
{
    if (record.cpu.time) {
        time.1 <- as.numeric(proc.time()[3])
    }
    else {
        time.1 <- "not timing"
    }
    if (method == "LM") {
        method.name = "Method of L-Moments"
    }
    if (method == "LM") {
        mx.results <- fit.gpd.lmom(data = x, na.rm = na.rm)
    }
    if (record.cpu.time) {
        time.2 <- as.numeric(proc.time()[3])
        runtime <- round(time.2 - time.1, 2)
    }
    else {
        runtime = NA
    }
    result <- list(lambda = mx.results[, 1], mx = mx.results, 
        cpu.time = runtime, param = "gpd", method.code = method, 
        method.name = method.name)
    if (return.data) {
        result$data = x
    }
    class(result) <- "starship"
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
