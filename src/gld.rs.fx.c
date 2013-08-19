/* gld.rs.fx.c - Part of Robert King's gld package for the R statistical
 * language.
 *
 * Copyright (C) Robert King 1993,2000,2001
 * robert.king@newcastle.edu.au
 * http://maths.newcastle.edu.au/~rking/publ/software.html
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
 * These functions calculate F(x) for the rs paramterisation for the gld
 * This is not very good C.  I wrote this in 1993, as my first real program
 * and I haven't tidied it up very much for release with the R package.
 * However, it has been used for some time and works fine.
 *
 * This is code is loaded into R.  It doesn't need stdio.h 
 * It needs math.h, which is loaded automagically by R.
 */

#include <stdio.h> 
#include <math.h> 
#include <R.h> 

static double la, lb, lc, ld, x;

/* the function that finds the root */

void gl_rs_distfunc( double *pa,double *pb,double *pc,double *pd, 
double  *px1,double *px2,double *pxacc, int *max_it,
double *ecks, double *u, int *pl)
{

/* pa to pd: 	pointers to the values of the parameters of the gld (rs param)
 * px1:		minimum value of u, should be zero
 * px2:		maximum value of u, should be 1
 * pxacc:	desired accuracy of the calculation
 * max_it: 	maximum iterations for N-R root finder 
 * ecks:	the quantiles of the gld given
 * u:		array to put the calculated depths
 * pl:		length of the data
 */
		
	double  x1, x2, xacc; 		
	double  a, b, c, d;
	int l;

	int i,j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	void funcd();

	x1 = *px1; x2 = *px2; xacc = *pxacc;
	a = *pa; b = *pb; c = *pc; d = *pd;
	l = *pl;

	la = a; lb = b; lc = c; ld = d;

/* The C version has something here to force the limits to be xacc and 1-xacc
	rather than 0 and 1 if lambda3 and lambda4 are negative.  I can't 
	see why, so I'm leaving it out.*/

for (i=0;i<l;i++)
{
        x = ecks[i];
	u[i] = 0.0;
	funcd(x1,&fl,&df);
	funcd(x2,&fh,&df);
	if (fl*fh >= 0.0) {
		error("gld package C code numerical failure (this should not happen - please report to maintainer)\n Program aborted during calculation of F(x)\n at parameter values %f, %f, %f, %f\n The x value was index: %d, value %f\n",*pa, *pb, *pc, *pd, i, x);
		}
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
		}
	else {
		xh = x1;
		xl = x2;
		}
	rts = 0.5*(x1+x2);
	dxold = fabs(x2-x1);
	dx = dxold;
	funcd(rts,&f,&df);
	for (j=1;j<= *max_it;j++) {
		if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0.0 ) ||
( fabs(2.0*f) > fabs (dxold*df))) {
			dxold = dx;
			dx = .5* (xh - xl);
			rts = xl +dx;
			if (xl == rts ) { u[i] = rts; break; }
			}
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) { u[i] = rts; break; }
			}
		if (fabs(dx) < xacc) { u[i] = rts; break; }
		funcd(rts,&f,&df);
		if (f < 0.0)
			xl =rts;
		else 
			xh =rts;
		}
}
}


void funcd(u,a,b)
double u,*a,*b;
{

/* function is the gld F-1(u)  */
/* df_du is dF-1(u)/du  */

double function, df_du;

if ( lc == 0 ) {
        if ( ld == 0 ) {
                /*      Both l3 and l4 zero     */
                function = la - x;
                df_du = 0;
                }
        else {
                /*      l3 zero, l4 non-zero    */
                function = la + ( 1 - pow((1-u),ld) )/ lb - x;
                df_du = (ld) * ( pow((1.0-u),(ld - 1.0)) / lb);
                }
        }
else {
        if ( ld == 0 ) {
                /*  l4 zero, l3 non-zero    */
                function = la + ( pow(u,lc) - 1 ) / lb - x;
                df_du = lc * ( pow(u,(lc- 1.0)) ) / lb;
                }
        else {
                /*  l3 non-zero, l4 non-zero    */
                function = ( pow((u),lc) - pow((1.0-u),ld) )/ lb + la - x;
                df_du = (lc * (pow((u),(lc-1.0))) + ld * ( pow((1.0-u),(ld-1.0)) ))/ lb;
                }
        }


*a = function;
*b = df_du;

}

