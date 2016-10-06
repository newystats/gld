# (c) 2016 Robert King.  GPL (>= 2)
# This function calculates the Beta Function of the same lambda.
# It is needed for the standard error calculations for the gld gpd.
# The gld package doesn't export it.
# This function is more than needed, really because the SE only exists
# for lambda > -0.5
BetaLambdaLambda <- function(lambda){
  some.dodgy <- FALSE # identifies if there are problem lambda values
  if (any(lambda[lambda<1]%%1==0)) { # If any of the lambdas are negative integers
    some.dodgy <- TRUE
  }
  if (any(lambda == -0.5)) { some.dodgy <- TRUE} # Beta(-0.5,-0.5) = 0, but you can't get it directly from the version involving the Gamma Function
  if (some.dodgy){
    dodgy.ones <- (((lambda <= 0) & (lambda%%1==0) ) | lambda== -0.5)
    ok.ones <- lambda[!dodgy.ones] # Finds the values of lambda that work for this version
    ok.result <- (2 * (gamma(ok.ones+1))^2) / (ok.ones *gamma(2*ok.ones + 1))
    result <- c()
    result[dodgy.ones] <- NaN # These could possibly be Inf, but the function loops round at negative integers so NaN seems a better option
    result[lambda == -0.5] <- 0
    result[!dodgy.ones] <- ok.result
  } else {
    result <- (2 * (gamma(lambda+1))^2) / (lambda *gamma(2*lambda + 1))
    # Follows from Beta fn as ratio of gammas and $$\Gamma(-1+\epsilon) = \frac{\Gamma(\epsilon) }{-1+\epsilon}$$
  }
  result
}