qgl <- function(p,lambda1,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL)
{
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the values are OK
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),"\ndo not produce a proper distribution for the",param,"parameterisation \n - see documentation for gl.check.lambda"))
	}
result <- switch(param,
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	FKML=,
	fkml=,
	fmkl=.qgl.fmkl(p,lambdas),
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs=.qgl.rs(p,lambdas),
	fm5 = .qgl.fm5(p, lambdas),
	stop("Error: Parameterisation must be fkml, fm5 or rs")
	) # closes "switch"
result
}


.qgl.fmkl <- function(p,lambdas)
{
# No checking - use qgl if you want that
lambda4 = lambdas[4]
lambda3 = lambdas[3]
lambda2 = lambdas[2]
lambda1 = lambdas[1]
p <- as.double(p)
# abandoned this for the simpler one below
# outside.range <- !as.logical(((p<1)*(p>0))|(sapply(p, all.equal,1)=="TRUE")| (sapply(p, all.equal, 0)=="TRUE"))
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# If OK, determine special cases
if (lambda3 == 0) { 
	if (lambda4 == 0) { # both log
		quants <- lambda1 + (log(u) - log(1 - u))/lambda2
		}
	else	{ # l3 zero, l4 non-zero
		quants <- lambda1 + 
			(log(u) - ((1 - u)^lambda4 - 1)/lambda4)/lambda2
		}
	}
else 	{ # lambda3 non-zero
	if (lambda4 == 0) { # non-zero, l4 zero
		quants <- lambda1 + 
			((u^lambda3 - 1)/lambda3 - log(1 - u))/lambda2
		}
	else	{ # both non-zero - use usual formula
		quants <- lambda1 + ( ( u ^ lambda3 - 1)	/ lambda3 
			- ( (1-u)^lambda4 - 1) /lambda4 ) / lambda2
		}
	}
# Now we have the quantiles for p values inside [0,1], put them in the right
# place in the result vector
result <- p*NaN
result[!outside.range] <- quants
# The remaining "quantiles" are NaN
result
}

.qgl.fm5 <- function(p,lambdas)
{
# No parameter value checking. If you want that, use qgl!
lambda5 = as.double(lambdas[5])
lambda4 = as.double(lambdas[4])             
lambda3 = as.double(lambdas[3])
lambda2 = as.double(lambdas[2])
lambda1 = as.double(lambdas[1])
p <- as.double(p)
# abandoned this for the simpler
# outside.range <- !as.logical(((p<1)*(p>0))|(sapply(p, all.equal,1)=="TRUE")| (sapply(p, all.equal, 0)=="TRUE"))
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# If OK, determine special cases
if (lambda3 == 0) { 
	if (lambda4 == 0) { # both log
		quants <- lambda1 + ((1-lambda5)*log(u) - (1+lambda5)*log(1 - u))/lambda2
		}
	else	{ # l3 zero, l4 non-zero
		quants <- lambda1 + 
			((1-lambda5)*log(u) - (1+lambda5)*((1 - u)^lambda4 - 1)/lambda4)/lambda2
		}
	}
else 	{ # lambda3 non-zero
	if (lambda4 == 0) { # non-zero, l4 zero
		quants <- lambda1 + 
			((1-lambda5)*(u^lambda3 - 1)/lambda3 - (1+lambda5)*log(1 - u))/lambda2
		}
	else	{ # both non-zero - use usual formula
		quants <- lambda1 + ((1-lambda5)* ( u ^ lambda3 - 1) / lambda3 
			- (1+lambda5)*( (1-u)^lambda4 - 1) /lambda4 ) / lambda2
		}
	}
# Now we have the quantiles for p values inside [0,1], put them in the right
# place in the result vector
result <- p*NaN
result[!outside.range] <- quants
# The remaining "quantiles" are NaN
result
}

.qgl.rs <- function(p,lambdas)
{
u <- p
# No parameter value checking - use qgl!
lambda4 = lambdas[4]
lambda3 = lambdas[3]
lambda2 = lambdas[2]
lambda1 = lambdas[1]
# At present, I'm rejecting zero values for l3 and l4, though I think there 
# are limit results, so one functional form.
quants <- lambda1 + ( u ^ lambda3 - (1-u)^lambda4 ) / lambda2
quants
}


qdgl <- function(p,lambda1,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL){
	dqgl(p=p,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3,lambda4=lambda4,param=param,lambda5=lambda5)}

dqgl <- function(p,lambda1,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL)
{
# Don't allow characters in lambda5 - common error with parameterisation stuff
if(is.character(lambda5)) {stop(paste("lambda5=",lambda5,"It should be a number between -1 and 1"))}
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the values are OK
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),"\ndo not produce a proper distribution for the",param,"parameterisation \n - see documentation for gl.check.lambda"))
	}
result <- switch(param,  
# Different tests apply for each parameterisation
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	FKML=,
	fkml=,
	fmkl=.dqgl.fmkl(p,lambdas),
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs=.dqgl.rs(p,lambdas),
	fm5 = .dqgl.fm5(p, lambdas),
	stop("Error: Parameterisation must be fkml, fm5 or rs")
	) # closes "switch"
result
}


.dqgl.rs <- function(p,lambdas)
{
# Check the values are OK)
if(!gl.check.lambda(lambdas,param="rs",vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),
"\ndo not produce a proper distribution with the RS parameterisation - see \ndocumentation for gl.check.lambda"))
	}
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]	
dens <-  lambdas[2]/(lambdas[3] * (p^(lambdas[3] -1)) + lambdas[4] * ((1 - p)^(lambdas[4] -1)))
dens
}


.dqgl.fmkl <- function(p,lambdas)
{
# Check the values are OK)
if(!gl.check.lambda(lambdas,param="fkml",vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),
"\ndo not produce a proper distribution with the FMKL",
"parameterisation - see \ndocumentation for gl.check.lambda"))
	}
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# The density is given by 1/Q'(u)
dens <- lambdas[2]/(p^(lambdas[3] - 1) + (1 - p)^(lambdas[4] - 1))
dens
}

.dqgl.fm5 <- function(p,lambdas)
{
# Check the values are OK)
if(!gl.check.lambda(lambdas,param="fm5",vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),
"\ndo not produce a proper distribution with the FM5",
"parameterisation - see \ndocumentation for gl.check.lambda"))
	}
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# The density is given by 1/Q'(u)
dens <- lambdas[2]/((1-lambdas[5])*(u^(lambdas[3] - 1)) + (1+lambdas[5])*((1 - u)^(lambdas[4] - 1)) )
dens
}

