# omnibus gpd fitter - modelled on fit.fkml
fit.gpd <- function(x,method="LM",na.rm=TRUE,
      record.cpu.time = TRUE, return.data=FALSE){
  # rec.cpu.time TRUE for consistency with fit.fkml
  
  # Start timing
  if(record.cpu.time) { 
    time.1 <- as.numeric(proc.time()[3]) } else { time.1 <- "not timing"
    }
  
  # set long method name
  if (method == "LM")  {method.name="Method of L-Moments"}
  
  if (method == "LM")  {
    mx.results <- fit.gpd.lmom(data=x,na.rm=na.rm)
  }
  # results will always be a matrix (with SEs)
  
  # Store results
  if (record.cpu.time) {time.2 <- as.numeric(proc.time()[3]); runtime <- round(time.2-time.1,2) } else {runtime=NA}
  
  result <- list(mx=mx.results)
  #names(result$lambda) <- paste("lambda",1:length(result$lambda),sep="")
  
   #result <- list(lambda= # a vector 
   #                unclass(pars.est[1:4]),grid.results=grid.results,
   #              optim.results=optim.results,param="fkml",
  #               method.code=method,method.name=method.name,
  #               fkml.oldstyle.results=old.style.results)
  if (return.data) {result$data = x}
  class(result) <- "starship"
   
  result <- "function not finished yet"
}
                   

fit.gpd.lmom <- function(data,na.rm=TRUE){
  if (na.rm){ dataNArm <- data[!is.na(data)] 
  } else { if (any(is.na(data))) {
      stop(paste("NA values in ",deparse(substitute(data)),". use na.rm=TRUE to fit these data.",sep=""))} else {dataNArm <- data}
  }
  fit.gpd.lmom.given(lmoms=lmom::samlmu(dataNArm,nmom=4),n=length(dataNArm))
}


fit.gpd.lmom.given <- function(lmoms,n=NULL){
  warning("calculations here are still sld")
  if (lmoms[3]>(1/3)) {stop("No QB skew logistic distribution corresponds to these L Moment values.\nThese L Moments are more right skew than the exponential distribution, the limiting case of the QB Skew Logistic.")} 
  if (lmoms[3]<(-1/3)) {stop("No QB skew logistic distribution corresponds to these L Moment values.\nThese L Moments are more left skew than the reflected exponential distribution, the limiting case of the QB Skew Logistic.")}
  ah = lmoms[1] - 6*(lmoms[3]*lmoms[2]) # alpha hat 
  bh = 2*lmoms[2] # beta hat
  dh = 0.5*(1+3*lmoms[3]) # delta hat
  lmomest <- c(ah,bh,dh)
  names(lmomest) <- c("alpha","beta","delta")
  if (!is.null(n)){ # Sample size is known - calculate Std Errors
    om = dh*(1-dh) # omega
    se.alpha = bh * sqrt((57 + (125*pi^2-1308)*om)/(15*n))
    se.beta = bh * sqrt(4/(3*n) * (1 - (pi^2-8)*om))
    se.delta = sqrt((8-(397+160*om-20*pi^2*(om+2))*om)/(15*n))
    lmomse <- c(se.alpha,se.beta,se.delta)
    ret <- cbind(lmomest,lmomse)
    dimnames(ret) <- list(c("alpha","beta","delta"),
                          c("Estimate","Std. Error"))
  } else {
    ret <- lmomest # return just the estimates
  }
  ret
}
