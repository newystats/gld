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
fit.gpd.lmom(data,na.rm=TRUE)
fit.gpd.lmom.given(lmoms,n=NULL)
}
\arguments{
  \item{x}{Data to be fitted, as a vector}
  \item{method}{A character string, to select the estimation method.  Only Method of
  L-Moments "LM" is implemented.}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{record.cpu.time}{Boolean --- should the CPU time used in fitting be recorded in the fitted model object?}
  \item{return.data}{Logical: Should the function return the data (from the argument \code{x})?}
}
\details{The method of L-Moments equates sample L-Moments with expressions for the 
L-Moments of the GPD type GLD.  Closed form expressions exist to give these estimates.

If 

For many values 
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
