# Standard errors of the gld gpd L moment estimates.
# These need the BetaLambdaLambda function that gives the special
# case of the beta extended to negative non integer lambda
# They also need the nu functions from gldgpdNus.R

se.alphahat <- function(alpha,beta,delta,lambda,n){
  omega <- delta*(1-delta)
  beta*sqrt(1/(n*(lambda-1)^2*(lambda+2)*nu1(lambda)))*sqrt(
               (nu2(lambda) - (omega/(lambda*(lambda+1)))*(
    (nu3(lambda)/4 - omega*nu4(lambda)*nu5(lambda))*BetaLambdaLambda(lambda)))+
      nu6(lambda)*(nu7(lambda)-4*omega*nu4(lambda)*nu8(lambda)))
}

se.betahat <- function(alpha,beta,delta,lambda,n){
  omega <- delta*(1-delta)
  beta*sqrt(1/(n*(lambda+2)*nu1(lambda))*(nu9(lambda)-(omega/(lambda*(lambda+1)))*(nu10(lambda)/4*BetaLambdaLambda(lambda)+nu6(lambda)*nu11(lambda))))
}

se.deltahat <- function(alpha,beta,delta,lambda,n){
  omega <- delta*(1-delta)
  sqrt(((lambda+1)*(lambda+2)*(lambda+3))/(n*(lambda-1)^2*nu1(lambda)))*(
    sqrt(nu12(lambda)-omega/lambda*(1/4*(nu13(lambda)-
    12*omega*nu14(lambda))*BetaLambdaLambda(lambda) - 
    nu6(lambda)/(lambda+3)*(3*nu15(lambda) + omega*nu16(lambda))) ))
}

se.lambdahat <- function(alpha,beta,delta,lambda,n){
  omega <- delta*(1-delta)
  sqrt(((lambda+1)*(lambda+2)*(lambda+3)^2*(lambda+4)^2)/(n*nu1(lambda))*(
    nu17(lambda) - omega/lambda*(1/4*nu5(lambda)*BetaLambdaLambda(lambda)+nu6(lambda)*nu8(lambda))))
}