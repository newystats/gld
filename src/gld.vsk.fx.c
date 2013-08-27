/* gld.vsk.fx.c - Part of Robert King's gld package for the R statistical
 * language.
 *
 * Copyright (C) Robert King 1993,2000,2001,2006,2011,2013
 * robert.king@newcastle.edu.au
 * http://tolstoy.newcastle.edu.au/~rking/publ/software.html
 *
 * This program is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) 
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA  02111-1307  USA
 *
 * These functions calculate F(x) for the van Staden and King parameterisation
 * of the gld.  It is adapted from gld.rs.fx.c (also in this R package).
 * 
 * This is code is loaded into R.
 * The three include lines are handled by R, but need to be here to keep the compile step happy
 */

#include <stdio.h>
#include <math.h>
#include <R.h> 

/* proto */

void vsk_funcd( double , double , double *, double *, double *, double *, double *, double *);

/* the function that finds the root: this is almost generic */
/* Since it is so generic, this function retains the parameters pa to pd, 
  rather than renaming for the vsk */

void gl_vsk_distfunc( double *pa,double *pb,double *pc,double *pd, 
double  *pu1,double *pu2,double *pxacc, int *max_it,
double *ecks, double *u, int *lengthofdata)
{

/* pa to pd:    pointers to the values of the parameters of the gld (vsk param)
 * pa		alpha
 * pb		beta
 * pc		delta
 * pd		lambda
 * pu1:         minimum value of u, usually zero
 * pu2:         maximum value of u, usually 1
 * pxacc:       desired accuracy of the calculation
 * max_it:      maximum iterations for N-R root finder
 * ecks:        the quantiles of the gld given
 * u:           array to put the calculated depths
 * pl:          length of the data
 */

double  u1, u2, xacc; 

int i,j; /* assorted indices */
double df,dx,dxold,f,fh,fl; /* to invert the equation */
double x;
double temp,xh,xl,rts;
/* initialising things */
i=0; j=0;
df=0.0;dx=0.0;dxold=0.0;f=0.0;fh=0.0;fl=0.0;
x=0.0;
temp=0.0;xh=0.0;xl=0.0;rts=0.0;

u1 = *pu1; u2 = *pu2; xacc = *pxacc;

if (*pc < 0) {  /* shouldn't happen for this parameterisation */
	if (u1 == 0) {
		u1 = xacc;
		}
	if (u2 == 1) {
		u2 = 1-xacc;
		}	
	}
if (*pd < 0) {
	if (u1 == 0) {
		u1 = xacc;
		}
	if (u2 == 1) {
		u2 = 1-xacc;
		}	
	}


for (i=0;i<*lengthofdata;i++)
{
    x = ecks[i];
	u[i] = 0.0;
	vsk_funcd(u1,x,&fl,&df,pa,pb,pc,pd);
	vsk_funcd(u2,x,&fh,&df,pa,pb,pc,pd);
	if (fl*fh >= 0.0) 
	{
		/* This is suggested in writing R extensions, but still gives the warning */
		error("Program aborted at parameter values %f, %f, %f, %f\n The data value being investigated was index %d, value: %f\n iteration %d, fl: %f fh: %f u1: %f u2: %f df: %f", *pa, *pb, *pc, *pd, i, x,j,fl,fh,u1,u2,df);
	}
	if (fl < 0.0) {
		xl = u1;
		xh = u2;
		}
	else {
		xh = u1;
		xl = u2;
		}
	rts = 0.5*(u1+u2);
	dxold = fabs(u2-u1);
	dx = dxold;
	vsk_funcd(rts,x,&f,&df,pa,pb,pc,pd);
	for (j=1;j<=*max_it;j++) {
		if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0.0 ) ||
( fabs(2.0*f) > fabs (dxold*df))) {
			dxold = dx;
			dx = .5* (xh - xl);
			rts = xl +dx;
			if (xl == rts ) { 
				u[i] = rts; 
				break; }
			}
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) { 
				u[i] = rts; 
				break; }
			}
		if (fabs(dx) < xacc) { 
			u[i] = rts; 
			break; }
		vsk_funcd(rts,x,&f,&df,pa,pb,pc,pd);
		if (f < 0.0)
			xl =rts;
		else 
			xh =rts;
		}
}
}


    void vsk_funcd( double u, double x, double *F, double *dFdu, double *alpha, double *beta, double *delta, double *lambda)
{
    /* double Qu,shouldnotbe0,shouldbe0;  for debugging */
/* *F is the gld F-1(u), with x subtracted */            
/* *dFdu is dF-1(u)/du	*/
/* NOTE: *dFdu is 1/f(x) a.k.a. 1/f(F-1(u)) - the opposite of what's required for the pdf */

if ( *lambda == 0 ) {
	/* lambda zero - special case */
    *F = *alpha + *beta * ( (1.0-*delta)*log(u) - *delta * log(1.0-u) ) - x;
    *dFdu = *beta * ( (1.0-*delta)/u + *delta /(1.0-u)) ;
    }
else {
    /*	lambda non-zero	*/
    /* Debugging
     shouldnotbe0 = (1.0-*delta)*(pow((u),*lambda) - 1.0 )/ *lambda;
    shouldbe0 = (*delta * (pow((1.0-u),*lambda) - 1.0))/ *lambda;
    Qu = *alpha + *beta * ( (1.0-*delta)*(pow((u),*lambda) - 1.0 )/ *lambda  - (*delta * pow((1.0-u),*lambda) - 1.0 )/ *lambda ); */
    *F = *alpha + *beta * ( (1.0-*delta)*(pow((u),*lambda) - 1.0 )/ *lambda  - (*delta * (pow((1.0-u),*lambda) - 1.0) )/ *lambda ) - x;
    *dFdu = *beta *( (1.0 - *delta)*( pow((u),(*lambda-1.0)) ) + ( *delta * pow((1.0-u),  (*lambda-1.0)) ) );
	}

}

