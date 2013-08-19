.gl.parameter.tidy <- function(lambda1,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL) 
{
# Don't allow characters in lambda5 - common error with parameterisation stuff
if(is.character(lambda5)) {stop(paste("lambda5=",lambda5,"It should be a number between -1 and 1"))}
# Don't allow numbers in parameterisation - included as a warning here, so the main one is a stop.
if(!is.character(param)) {warning(paste("param=",param,"It shouldn't be a number, it should be a string describing the parameterisation"))}
if(is.null(lambda1)) { stop("No values provided for lambda parameters, argument lambda1 is NULL") }
if(length(lambda1) > 1) #using a vector for the parameters.  
	# Check that there aren't values in the individual lambda arguments
	{
	if (!(is.null(lambda2) & is.null(lambda3)& is.null(lambda4) & is.null(lambda5)) ) 
		{ stop("Call includes vector version of the lambda parameters as well as the \nscalar version") }
	if ((length(lambda1) < 4) | (length(lambda1) > 5 ) )  
		{ stop(paste("argument lambda1 has length", length(lambda1),"\nThis should be 1 (lambda parameters given as seperate arguments), 4 (vector argument \n for RS or FKML parameterisation) or 5 (vector argument for fm5 parameterisation")) }
	if (length(lambda1)== 5)
		{ if (param != "fm5") { 
			stop(paste("argument lambda1 has length",length(lambda1),"which is not valid for the",param,"\nparameterisation")) 
			}
		# else --- fm5, in vector form, ready for gl.check.lambda 
		}
	if (length(lambda1)== 4)
		{ if (param == "fm5" ) 
			{ stop(paste("argument lambda1 has length 4, which is not valid for the fm5 \nparameterisation")) }
		# else --- 4 parameter versions in vector form, ready for gl.check.lambda 
		}
	}
else { # single parameter arguments - check they are there, then collect them together
	if (is.null(lambda2)) { stop("No value for lambda2") }
	if (is.null(lambda3)) { stop("No value for lambda3") }
	if (is.null(lambda4)) { stop("No value for lambda4") }
	if ((is.null(lambda5)) & param=="fm5" ) { stop("No value for lambda5") }
	if (!(is.null(lambda5)) & param!="fm5") { stop(paste("lambda5=",lambda5," but there is no lambda 5 for the\n",param,"parameterisation")) }
	if (param != "fm5") { # A 4 parameter version
		lambda1 <- c(lambda1,lambda2,lambda3,lambda4)
		}
	else { # fm5
		lambda1 <- c(lambda1,lambda2,lambda3,lambda4,lambda5)
		}
	}
# There is now an error if there is the wrong number of parameters, and 
# lambda1 returned as a vector with 4 or 5 elements
# as.double is needed to remove data.frame attributes if lambda1 was
# extracted from a data.frame
as.double(lambda1)
}

gl.check.lambda <-  function(lambdas,lambda2=NULL,lambda3=NULL,lambda4=NULL,
param="fkml",lambda5=NULL,vect=FALSE)
{
# Checks to see that the lambda values given are allowed.
# There is a function called .gl.parameter.tidy that does the tidying 
# around of parameters.  It return a single vector, which contains the
# parameters.
# If you call this after .gl.parameter.tidy, let it know with the vect=T
# argument
# If vect=TRUE, we don't need to tidy
if (vect) {
	if (!is.null(lambda3)) {
	warning("lambda3 should be null because you claim the parameters are in a vector")
		}
	}
else	{
	lambdas <- .gl.parameter.tidy(lambdas,lambda2,lambda3,lambda4,param,lambda5)
	}
if(param=="fm5"){lambda5 = lambdas[5]}
lambda4 = lambdas[4]                        
lambda3 = lambdas[3]
lambda2 = lambdas[2]
lambda1 = lambdas[1]
# I did have a check for finite lambdas, but that caused a problem with data frames, 
# so I removed it - still need to include the limit results
param <- switch(param,  
# Different tests apply for each parameterisation
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	FKML=,
	fmkl=,
	fkml={
	if (lambda2<=0) {return(FALSE)}
	else {return(TRUE)}
	},
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs={
	if (lambda3*lambda4>0) { # regions 3 and 4 
				 # all values of lambda 3 and lambda 4 OK
				 # check lambda 2
		if ((lambda3>0)&(lambda4>0)) { # region 3 - l2 >0
			if (lambda2<=0) {ret <- FALSE}
			else {ret <- TRUE}
			}
		if ((lambda3<0)&(lambda4<0)) { # region 4 - l2 <0
			if (lambda2>=0) {ret <- FALSE}
			else {ret <- TRUE}
			}
		}	
	else { 	# other quadrants - lambda 2 must be negative, and lambda3 
		# lambda 4 values need checking.
		if (lambda2>=0) {return(FALSE)}
		# Rectangular regions where RS is not defined 
		if ((lambda3>0)&(lambda3<1)&(lambda4<0)) {return(FALSE)}
		if ((lambda4>0)&(lambda4<1)&(lambda3<0)) {return(FALSE)}
		# Different here because there are a 
		# number of ways in which the parameters can fail.
		# 
		# Curved regions where RS is not defined
		# change to shorter var names
		lc <- lambda3
		ld <- lambda4
		if ((lambda3>-1)&(lambda3<0)&(lambda4>1)) {  # region 5 or not?
			if ( ((1-lc)^(1-lc)*(ld-1)^(ld-1))/((ld-lc)^(ld-lc)) > -lc/ld )	
				{return(FALSE)}
			else 	{return(TRUE)}
			}
		# Second curved region 
		if ((lambda4>-1)&(lambda4<0)&(lambda3>1)) {  # region 6 or not?
			if ( ((1-ld)^(1-ld)*(lc-1)^(lc-1))/((lc-ld)^(lc-ld)) > -ld/lc )
				{return(FALSE)}
			else 	{return(TRUE)}
			}
		# There may be some limit results that mean these are not correct, but
		# I'll check that later
		# This is not the place where the possible l3,l4 zero values should appear
		if (lambda3 == 0) {
			warning('lambda 3 = 0 with RS parameterisation - possible problem')
			if (lambda4 == 0) {return(FALSE)}
			else {return(TRUE)}
			}
		if (lambda4 == 0) {
			warning('lambda 5 = 0 with RS parameterisation - possible problem')
			if (lambda4 == 0) {return(FALSE)}
			else {return(TRUE)}
			}
		# If we get here, then the parameters are OK.
		ret <- TRUE
		}
	},
	fm5={
		# make lambda5 - in here so it doesn't stuff up the other parameterisations
		lambda5 <- lambdas[5]
		if (lambda2<=0) {ret <- FALSE}
		else { #Check lambda5 - should be between -1 and 1, but I haven't checked this against a piece of paper
			if ((lambda5 >= -1) | (lambda5 <= 1)) {ret <- TRUE}
			else {ret <- FALSE}
		}	
	},
	stop("Error when checking validity of parameters.\n Parameterisation must be fmkl, rs or fm5")
	) # closes "switch"
ret
} 
