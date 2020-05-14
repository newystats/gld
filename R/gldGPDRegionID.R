lambdaTilde = sqrt(6)-1 # boundary between regions A & B. Not exported

gldGPDRegionID = function(pars=NULL,lambda=NULL){ # identify the region. Not exported 
  if (is.null(lambda)) {
    if (is.null(pars)) {stop("no values given for either pars or lambda")
    } else {
      lambda = pars[4]
    }
  } else {
    if (is.null(pars)) {
      # all good
    } else { #check that lambda is the same as the 
      if (lambda != pars[4]) {
        warning(paste("both pars (",paste(pars,collapse = ", "),") and lambda (",lambda,") given with different lambda \nvalues. Using value from pars:",pars[4],collapse=""))
        lambda = pars[4]
      }
    }
  }
  if (lambda < lambdaTilde){
    region = "A"
  } else {
    region = "B"
  }
  region
}
