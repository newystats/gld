# gld
R package for the Generalised Lambda Distribution.  This is available on CRAN as https://cran.r-project.org/package=gld

The generalised $\lambda$ distribution, also known as the Asymmetric Lambda or Asymmetric Tukey Lambda Distribution provides a wide variety of shapes with one functional form.   

This package provides a number of different _types_ of the distribution, each defined by their [Quantile function](https://en.wikipedia.org/wiki/Quantile_function):

1. FMKL type (default) -

$$Q(u)= \lambda_1 + { { \frac{u^{\lambda_3}-1}{\lambda_3} - 	
\frac{(1-u)^{\lambda_4}-1}{\lambda_4} } \over \lambda_2 } 
\qquad\text{ for } \quad\lambda_2 > 0$$

[Freimer Mudholkar, Kollia and Lin (1988)](https://doi.org/10.1080/03610928808829820) 
