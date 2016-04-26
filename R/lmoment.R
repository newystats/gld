gld.lmoments <- function(pars,order=1:4,ratios=TRUE,type="GPD"){
  # pars is the vector of parameters for the GPD type GLD
  gl.check.lambda(lambdas = pars,param = type)
  if (length(order) == 0) {stop ("argument order of gld.lmoments requires values")}
  if (any(is.na(order))) {stop ("argument order of gld.lmoments does not support NA values")}
  order.orig <- order
  order <- round(order.orig)
  if ( min(order) < 1) {warning("At least one element of the order argument is less than 1.  This function implements L Moments of positive whole number order only")}
  if ( max(abs(order - order.orig)) > .Machine$double.eps^0.5 ){
    warning("At least one element of the order is not an integer.  This function implements L Moments of positive whole number order only")
  }
  alpha = pars[1]
  beta = pars[2]
  delta = pars[3]
  lambda = pars[4]
  if (lambda <= -1) {return(rep(Inf,length(order)))}
  if (max(order)>4){
    # other than first 4
    # match 1,2 in order and do those, then use standard formulae
    L1 = alpha + (beta*(2*delta -1))/(lambda +1 )
    L2 = beta / ((lambda+1)*(lambda+2))
    result <- rep(NA,length(order))
    names(result) <- paste("L",order,sep="")
    result[(order==1)] <- L1
    result[(order==2)] <- L2
    if (ratios){
    result[(order>2)]<-((gamma(3+lambda))*(gamma(order[(order>2)]-1-lambda)))/
                        ((gamma(1-lambda))*(gamma(order[(order>2)]+1+lambda)))
    taunames <- paste("tau",order,sep="")
    names(result)[(order>2)] <- taunames[(order>2)]
    } else {
      result[(order>2)]<-((gamma(1+lambda))*(gamma(order[(order>2)]-1-lambda)))/
        ((gamma(1-lambda))*(gamma(order[(order>2)]+1+lambda)))
      }
  } else {
  L1 = alpha + (beta*(2*delta -1))/(lambda +1 )
  L2 = beta / ((lambda+1)*(lambda+2))
  T3 = ((2*delta -1)*(1-lambda))/(lambda+3)
  T4 = ((lambda-1)*(lambda-2))/((lambda+3)*(lambda+4))
  if (ratios){
    result <- c(L1,L2,T3,T4)
    names(result) <- c("L1","L2","tau3","tau4")
    result <- result[order]
  } else {
    result <- c(L1,L2,T3*L2,T4*L2)
    names(result) <- paste("L",1:4,sep="")
    result <- result[order]
    }
  }
  result
}