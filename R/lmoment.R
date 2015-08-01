LmomentGLD <- function(lambdas,order=1:4,ratios=TRUE,param="GPD"){
  # lambda is the vector of lambda parameters for the FKML type
  gl.check.lambda(lambdas = lambdas)
  order.orig <- order
  order <- round(order.orig)
  if ( min(order) < 1) {warning("At least one element of the order argument is less than 1.  This function implements L Moments of positive whole number order only")}
  if ( max(abs(order - order.orig)) > .Machine$double.eps^0.5 ){
    warning("At least one element of the order is not an integer.  This function implements L Moments of positive whole number order only")
  }
  alpha = lambdas[1]
  beta = lambdas[2]
  delta = lambdas[3]
  lambda = lambdas[4]
  if (lambda <= -1) {return(rep(Inf,length(order)))}
  L1 = alpha + beta*(2*delta -1)/(lambda +1 )

  L1
}