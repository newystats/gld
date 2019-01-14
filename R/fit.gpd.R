
# omnibus gpd fitter - modelled on fit.fkml
fit.gpd <- function(x,method="LM",na.rm=TRUE,
      record.cpu.time = FALSE, return.data=FALSE){
  # rec.cpu.time TRUE for consistency with fit.fkml
  
  # Start timing
  if(record.cpu.time) { 
    time.1 <- as.numeric(proc.time()[3]) } else { time.1 <- "not timing"
    }
  
  # set long method name
  if (method == "LM")  {method.name="Method of L-Moments"}
  
  if (method == "LM")  {
    results <- fit.gpd.lmom(data=x,na.rm=na.rm)
  }
  # results a list with 2 elements
  
  # add optional elements
  if (record.cpu.time) {
    time.2 <- as.numeric(proc.time()[3]); runtime <- round(time.2-time.1,2) 
    results$cpu <- runtime 
    }
  if (return.data) {results$data = x}
  class(results) <- "GldGPDFit"
  results
}
                   

fit.gpd.lmom <- function(data,na.rm=TRUE){
  if (na.rm){ dataNArm <- data[!is.na(data)] 
  } else { if (any(is.na(data))) {
      stop(paste("NA values in ",deparse(substitute(data)),". use na.rm=TRUE to fit these data.",sep=""))} else {dataNArm <- data}
  }
  fit.gpd.lmom.given(lmoms=lmom::samlmu(dataNArm,nmom=4),n=length(dataNArm))
}

fit.gpd.lmom.given <- function(lmoms,n=NULL){
  if (length(lmoms) < 4) {stop("4 L-Moments are required to fit the GLD gpd.\nArgument lmoms of fit.gpd.lmom.given is less than 4 long.")}
  t4 <- lmoms[4]
  t3 <- lmoms[3]
  l2 <- lmoms[2]
  l1 <- lmoms[1]
  el.1 <- (3+7*t4)
  if (abs(t3)>=1){problem=paste("No estimates possible, impossible sample Tau 3 value: Tau3=",t3,"outside (-1,1) range")
    warning(problem)
    res <- list(estA=NA,estB=NA,warn=problem,param="gpd")
    class(res) <- "GldGPDFit"
    return(res)}
  if ( (5*t3^2-1)/4 > t4 ){problem = paste("No estimates possible, impossible sample Tau3/Tau4 combination. (5*Tau3^2-1)/4 =",(5*t3^2-1)/4,"must be <= Tau4 =",t4)
    warning(problem)
    res <- list(estA=NA,estB=NA,warn=problem,param="gpd")
    class(res) <- "GldGPDFit"
    return(res)}
  if (t4>=1){problem = paste("No estimates possible, impossible sample Tau 4 value: Tau4=",t4,">= 1")
    warning(problem)
    res <- list(estA=NA,estB=NA,warn=problem,param="gpd")
    class(res) <- "GldGPDFit"
    return(res)}
  if ((t4^2+98*t4+1)<0) {problem = paste("No estimates possible, Tau4 too low (lowest possible value is approx -0.0102051). Tau4 here is ",t4)
    warning(problem)
    res <- list(estA=NA,estB=NA,warn=problem,param="gpd")
    class(res) <- "GldGPDFit"
    return(res)}
  el.2 <- sqrt(t4^2+98*t4+1)
  denom <- (2*(1-t4))
  lambdahatA <- (el.1 - el.2 )/ denom
  lambdahatB <- (el.1 + el.2 )/ denom
  deltahatA <- 0.5*(1-(t3*(lambdahatA+3))/(lambdahatA-1))
  deltahatB <- 0.5*(1-(t3*(lambdahatB+3))/(lambdahatB-1))  
  betahatA <- l2*(lambdahatA+1)*(lambdahatA+2)
  betahatB <- l2*(lambdahatB+1)*(lambdahatB+2)
  alphahatA <- l1+(betahatA*(1-2*deltahatA))/(lambdahatA+1)
  alphahatB <- l1+(betahatB*(1-2*deltahatB))/(lambdahatB+1)
  if (isTRUE(all.equal(lambdahatA,0))) { # lamdbahatA is close to zero
    # Use SLD special case
    AisSLD = TRUE
    lambdahatA = 0 
  } else {AisSLD = FALSE}
  lmomestA <- c(alphahatA,betahatA,deltahatA,lambdahatA)
  if (gl.check.lambda(lmomestA,param="gpd")) {  
    names(lmomestA) <- c("alpha","beta","delta","lambda")
    RegionAest = TRUE} else {lmomestA <- NA
    RegionAest = FALSE}
  lmomestB <- c(alphahatB,betahatB,deltahatB,lambdahatB)
  if (gl.check.lambda(lmomestB,param="gpd")) {  names(lmomestB) <- c("alpha","beta","delta","lambda")} else {lmomestB <- NA}
 if (FALSE){ # code which is iffed out below
   # calculates SEs in the SLD special case, but not in 
   # the main gld GPD case if (!is.null(n)){ # Sample size is known - calculate Std Errors
    # Calculate standard errors for gld
    # Check if SEs exist in region A (lambda > -0.5 )
    #om = dh*(1-dh) # omega
    #se.alpha = bh * sqrt((57 + (125*pi^2-1308)*om)/(15*n))
    #se.beta = bh * sqrt(4/(3*n) * (1 - (pi^2-8)*om))
    #se.delta = sqrt((8-(397+160*om-20*pi^2*(om+2))*om)/(15*n))
    if (RegionAest){
      # calculate SEs for region A estimate
      # Don't do this test for region B, because they always exist?
      if (AisSLD){
        # lambda = 0, so this is the SLD - SE calculation ...
        # see fit.R in R package, sld
        warning("Since lambda estimate is zero, the estimated distribution\nis a special case, the Quantile Based Skew Logistic Distribution.\nNo standard errors are available for lambda, but SEs for the other\nparameters are given from the Quantile Based SLD.")
        omega = lmomestA[3]*(1-lmomestA[3])
        se.alpha = lmomestA[2] * sqrt((57 + (125*pi^2-1308)*omega)/(15*n))
        se.beta = lmomestA[2] * sqrt(4/(3*n) * (1 - (pi^2-8)*omega))
        se.delta = sqrt((8-(397+160*omega-20*pi^2*(omega+2))*omega)/(15*n))
        lmomestA <- cbind(lmomestA,c(se.alpha,se.beta,se.delta,NA))
        dimnames(lmomestA) <- list(c("alpha","beta","delta","lambda"),c("Estimate","Std. Error"))
        
      } else {warning("SE calculation not yet implemented")
      }
    } 
      # calculate SEs for region B estimate
      # If we have go this far, there is always a region B estimate?
      if (lmomestB[[4]] == 0){
        # lambda = 0, so this is the SLD - SE calculation ...
        # see fit.R in R package, sld
        warning("Since lambda estimate is zero, the estimate is a special case,\nthe Quantile Based Skew Logistic Distribution.  No standard errors are available for lambda,\nbut SEs for the other parameters are given from the Quantile Based SLD.")
        omega = lmomestB$delta*(1-lmomestB$delta)
        se.alpha = lmomestB$beta * sqrt((57 + (125*pi^2-1308)*omega)/(15*n))
        se.beta = lmomestB$beta * sqrt(4/(3*n) * (1 - (pi^2-8)*omega))
        se.delta = sqrt((8-(397+160*omega-20*pi^2*(omega+2))*omega)/(15*n))
        lmomestB <- cbind(lmomestB,c(se.alpha,se.beta,se.delta,NA))
        dimnames(lmomestB) <- list(c("alpha","beta","delta","lambda"),c("Estimate","Std. Error"))
        
      } else {warning("SE calculation not yet implemented")
      }
  } else {
    # return just the estimates
  }
  ret <- list(estA=lmomestA,estB=lmomestB,param="gpd")
  class(ret) <- "GldGPDFit"
  ret
}


