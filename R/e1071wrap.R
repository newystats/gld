# use e1071 to calculate moments for estimation
# this assumes na.rm has already happened
# don't export

moments4data <- function(x, type=1,var.divisor="n-1"){
  mean <- mean(x)
  if (var.divisor=="n"){
    n = length(x)
    var = ((n-1)/n)*var(x)
  } else { if (var.divisor=="n-1"){
    # n-1 
    var = var(x)
  } else {
    warning(paste("var.divisor should be n or n-1.  n-1 used instead of given value of",var.divisor))
    var = var(x)
    }
  }
  if (type==1){
  skewr <- e1071::moment(x,order = 3,center = TRUE, type=type)/(var^(1.5))
  kurtr <- e1071::moment(x,order = 4,center = TRUE, type=type)/(var^(2)) 
  # this version doesn't subtract 3
  } else {
    skew = e1071::skewness(x,type=type)
    kurt = e1071::kurtosis(x,type=type)+3 # note the +3, the gld fitting stuff uses the kurtosis, not the "excess kurtosis"
  }
  c(mean,var,skewr,kurtr)
}