# gld
R package for the Generalised Lambda Distribution.  This is available on CRAN as https://cran.r-project.org/package=gld

The generalised $\lambda$ distribution, also known as the Asymmetric Lambda or Asymmetric Tukey Lambda Distribution provides a wide variety of shapes with one functional form.   

This package provides a number of different _types_ of the distribution, each defined by their [Quantile function](https://en.wikipedia.org/wiki/Quantile_function):

1. FMKL type (default) -

$$Q(u)= \lambda_1 + { { \frac{u^{\lambda_3}-1}{\lambda_3} - 	
\frac{(1-u)^{\lambda_4}-1}{\lambda_4} } \over \lambda_2 } 
\qquad\text{ for } \quad\lambda_2 > 0$$

[Freimer Mudholkar, Kollia and Lin (1988) A study of the generalized tukey lambda family, Communications 
	in Statistics - Theory and Methods *17* 3547--3567.](https://doi.org/10.1080/03610928808829820) 

2. RS type -

$$ Q(u)= \lambda_1 + \frac{u^{\lambda_3} - (1-u)^{\lambda_4}}{\lambda_2} $$

[Ramberg, J. S. & Schmeiser, B. W. (1974) An approximate method for
generating asymmetric random variables, Communications of the ACM _17_ 78--82](https://dl.acm.org/doi/10.1145/360827.360840)

This type has a complex complex series of rules determining which values 
of the parameters produce valid statistical distributions.  See the gl.check.lambda function in this package for details.

3. GPD type -

$$ Q(u) = \alpha + \beta \left((1-\delta)\frac{(u^\lambda -1)}{\lambda} - \delta\frac{((1-u)^\lambda -1)}{\lambda} \right) 
\quad\text{ where }\beta > 0 \text{ and }-1 \leq \delta \leq 1$$

[Van Staden, Paul J., & M.T. Loots. (2009), Method of L-moment Estimation for the Generalized Lambda Distribution. In Proceedings of the Third Annual ASEARC Conference](https://eis.uow.edu.au/content/groups/public/@web/@inf/@math/documents/doc/uow074435.pdf)

4. FM5 type -

$$ Q(u) = \lambda_1 + \frac{ \frac{(1-\lambda_5)(u^{\lambda_3}-1)}{\lambda_3} - \frac{(1+\lambda_5)((1-u)^{\lambda_4}-1)}{\lambda_4} }{ \lambda_2 }$$
