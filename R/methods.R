print.starship <- function(x,digits = max(3, getOption("digits") - 3), ...)
{
# I should include the call
# Add names to the vector in starship
cat(paste(x$method.name,"estimate, gld type:",x$param,"\n"))
print.default(format(x$lambda,digits=digits), print.gap = 2,quote=FALSE)
}

summary.starship <- function(object,...)
{
cat(paste("Generalised Lambda Distribution",object$param,"type.",
          object$method.name," estimate.\n"))
cat("\nAdaptive Grid estimates:\n")
fake.lambda <- object$grid.results$lambda
names(fake.lambda) <- paste("lambda",1:length(fake.lambda),sep="")
fake.starship.object <- list(param=object$param,lambda=fake.lambda)
print.starship(fake.starship.object)
cat(paste("internal g-o-f measure at grid minimum:",
format(object$grid.results$response),"\n"))
cat("\nOptim (final) estimates (starting from grid estimates):\n")
print.starship(object)
cat(paste("internal g-o-f measure at optim minimum:",
format(object$optim.results$value),"\n"))
cat("optim.details:\nCounts: ")
print(object$optim.results$counts)
cat("Convergence: ")
print(object$optim.results$convergence)
cat("Message: ")
print(object$optim.results$message)
}

plot.starship <- function(x,data=NULL,ask=NULL,one.page=NULL,breaks="Sturges",plot.title="default",...)
{
if (plot.title == "default") {
  plot.title <- paste(x$method.name,"fit of",x$param,"type GLD")
}
allpar <- par()
opar <- allpar[match(c("ask","mfrow"),names(allpar))]
if (is.null(x$data)){
	if (is.null(data)) {stop("No data to compare to fit.  Use return.data=TRUE")} 
} else {
	if (is.null(data)) {data <- x$data #using data returned by starship function
		} else { 
		warning(paste(substitute(x),"has a data element and the data argument was also given.\nUsing ",paste(substitute(data))," instead of the data element of ",substitute(x))) } }
if (is.null(one.page)) { if (interactive()) {one.page=TRUE} else {one.page=FALSE}}
if (is.null(ask)) {if (one.page) {ask=FALSE} else {
  if (interactive()) {ask=TRUE} else {ask=FALSE}
    }
  }
if (ask) {par(ask=TRUE)}
if (one.page) {par(mfrow=c(2,1))}
qqgl(y=data,lambda.pars1=x$lambda,param=x$param,xlab=paste(x$method.name," Fitted Theoretical Quantiles"),main=plot.title) # add which option here
hist(data,prob=TRUE,xlab="Data",breaks=breaks,main=plot.title,...)
plotgld(lambda1=x$lambda,param=x$param,new.plot=FALSE,...)
par(opar) # Return to previous par
}

print.GldGPDFit <- function(x,digits = max(3, getOption("digits") - 3), ...)
{
  if (is.null(x$estA)) {
    if (is.null(x$estB)) {
      cat("No estimates for the GLD GPD\n")
      cat(x$warn)
    } else {
      cat("Region B only:\n")
      print.default(format(x$estB,digits=digits), print.gap = 2,quote=FALSE)
      # This needs to be extended once the package calculates SEs
    }
  } else { # region A estimate exists
    if (is.null(x$estB)){
      cat("Region A only:\n")
      print.default(format(x$estA,digits=digits), print.gap = 2,quote=FALSE)
    } else {
      cat("Region A:\n")
      print.default(format(x$estA,digits=digits), print.gap = 2,quote=FALSE)
      cat("\nRegion B:\n")
      print.default(format(x$estB,digits=digits), print.gap = 2,quote=FALSE)
    }
  }
}