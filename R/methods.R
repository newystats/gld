print.starship <- function(x,digits = max(3, getOption("digits") - 3), ...)
{
# I should include the call
# Add names to the vector in starship
cat(paste(x$method.name,"estimate, gld type:",x$param,"\n"))
if(!is.null(x$trim)){cat(paste("trimmed",x$trim[1],"low,",x$trim[2],"high. n=",x$trim[3],"\n"))}
print.default(format(x$lambda,digits=digits), print.gap = 2,quote=FALSE)
}

summary.starship <- function(object,...)
{
cat(paste("Generalised Lambda Distribution",object$param,"type.",
          object$method.name," estimate.\n"))
if(!is.null(object$trim)){
  cat(paste("L Moments based on data trimmed by;\n",object$trim[1],"lowest and",object$trim[2],"highest observations out of",object$trim[3]))
}
if(exists("object$grid.results")){
  cat("\nAdaptive Grid estimates:\n")
  fake.lambda <- object$grid.results$lambda
  if (object$method.code=="TL"){
      names(fake.lambda) <- paste("lambda",3:4)
    } else {
      names(fake.lambda) <- paste("lambda",1:length(fake.lambda),sep="")
    }
  fake.starship.object <- list(param=object$param,lambda=fake.lambda)
  print.starship(fake.starship.object)
  cat(paste("internal g-o-f measure at grid minimum:",
  format(object$grid.results$response),"\n"))
  cat("\nOptim (final) estimates (starting from grid estimates):\n")
} else {
  cat("\nOptim (final) estimates:\n")
}
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

plot.starship <- function(x,data=NULL,ask=NULL,one.page=FALSE,breaks="Sturges",plot.title="default",granularity=10000,...)
{
if (plot.title == "default") {
  plot.title <- paste(x$method.name,"fit of",x$param,"type GLD")
}
print("Using plot.starship from methods.R")
allpar <- par()
opar <- allpar[match(c("ask","mfrow"),names(allpar))]
if (is.null(x$data)){
	if (is.null(data)) {stop("No data to compare to fit.  Use return.data=TRUE")} 
} else {
	if (is.null(data)) {data <- x$data #using data returned by starship function
		} else { 
		warning(paste(substitute(x),"has a data element and the data argument was also given.\nUsing ",paste(substitute(data))," instead of the data element of ",substitute(x))) } }
# if (is.null(one.page)) { if (interactive()) {one.page=TRUE} else {one.page=FALSE}} # this was here so that interactive sessions didn't have a user-visible change
if (is.null(ask)) {if (one.page) {ask=FALSE} else {
  if (interactive()) {ask=TRUE} else {ask=FALSE}
    }
  }
if (ask) {par(ask=TRUE)}
if (one.page) {par(mfrow=c(2,1))}
qqgl(y=data,lambda.pars1=x$lambda,param=x$param,xlab=paste(x$method.name," Fitted Theoretical Quantiles"),main=plot.title) # add which option here
# do hist 2nd to fix problem with cutting off the top of the density? - doesnt 
# work - manually calculate max for density and adjust hist limits
#print("Which has this thing to make granularity.  Why are you complaining?")
default.granularity = 10000
#if(exists(granularity)) {
  if(is.integer(granularity)){
    if(granularity < 30){
      # too small
      granularity = default.granularity
    } else {
      if (granularity > 10000000) {
        # too big
        granularity = default.granularity
      }
    }
  } else {granularity = default.granularity}
u <- seq(from = 0, to = 1, by = 1/granularity)
# quantiles not neeeded quantiles <- qgl(u,lambda1=lambdas,param=param)
density <- dqgl(u,lambda1=x$lambda,param=x$param)
max.density = max(density)
hist(data,prob=TRUE,xlab="Data",breaks=breaks,ylim=c(0,max.density),...)  
plotgld(lambda1=x$lambda,param=x$param,main=plot.title,new.plot=FALSE,...)
if (one.page) {par(opar)} # Return to previous par
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

plot.GldGPDFit <- function(x,data=NULL,ask=NULL,breaks="Sturges",plot.title="default",col1="darkorange",col2="purple",...)
{
  if (plot.title == "default") {
    plot.title <- paste(x$method.name,"fit of",x$param,"type GLD")
  }
  # No one.page option here 
  allpar <- par()
  opar <- allpar[match(c("ask","mfrow"),names(allpar))]
  if (is.null(x$data)){
    if (is.null(data)) {stop("No data to compare to fit.  Use return.data=TRUE")} 
  } else {
    if (is.null(data)) {data <- x$data #using data returned by fit.gpd function
    } else { 
      warning(paste(substitute(x),"has a data element and the data argument was also given.\nUsing ",paste(substitute(data))," instead of the data element of ",substitute(x))) } }
  if (is.null(ask)) {
    if (interactive()) {ask=TRUE} else {ask=FALSE}
  }
  if (ask) {par(ask=TRUE)}  ## up to here ## 
  qqgl(y=data,lambda.pars1=x$lambda,param=x$param,xlab=paste(x$method.name," Fitted Theoretical Quantiles"),main=plot.title) # add which option here  
  # re-ordered to avoid chopping the top off the density - check this works
  hist(data,prob=TRUE,xlab="Data",breaks=breaks,main=plot.title,...)
  plotgld(lambda1=x$lambda,param=x$param,new.plot=FALSE,...)
  # No one.page option - check this - but I think we don't need to save par further up this function  if (one.page) {par(opar)} # Return to previous par
}
