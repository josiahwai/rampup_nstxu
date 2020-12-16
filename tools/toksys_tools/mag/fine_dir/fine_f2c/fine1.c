/*************************************************************************
SYNTAX: [hr,hz,fl] =fine1(ri,ro,zl,zu,cur,rp,zp,(idebug))
 
PURPOSE:  
   Magnetics function to determine the fields from an 
   axisymmetric rectangular uniform current region 
   (ri,ro,zl,zu,cur) to a point (rp,zp).
   routine operates on matrix of data: ri,ro,zl,zu,cur,rp,zp
 
INPUTS:
            ri,ro,zl,zu = inner,outer,lower,upper dimensions of coil [matrix]
            cur         = current [matrix]
            rp,zp       = receiving point [matrix]
            idebug      = debug parameter, 0=no, >0 = debug print out
OUTPUTS:
            hr, hz = fields at points rp,zp [Amp/m]
            fl     = flux at point [Vs/mu_o]
                     mu_o= 4*pi*1e-7
 

RESTRICTIONS:  This model is NOT VALID for currents approaching 0.
 
METHOD:  To compile in matlab: mex fine1.c fine.c -lm
 
WRITTEN BY:  Jim Leuer		ON 	4/?/95
Converted to C by: Mike Walker 	ON 	11/06/06
****************************************************************************/
#include "mex.h"
#include <math.h>

/* Input Arguments */

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
/*
     (all data is real*4 in fine but real*8 in matlab)

     pointers to MATLAB input real*8 arrays (PRHS)
*/

   double *ri_pt, *ro_pt, *zl_pt, *zu_pt, *cur_pt, *rp_pt, *zp_pt;
   double *idebug_pt;
   int k;

/*     debug parameter */
   int idebug;

/*
    size of the input matrix m=row, n=column, mn=total
*/
   int m, n, mn;

/*
    Left Hand Side Pointers (PLHS)
*/
   double *hr_pt, *hz_pt, *fl_pt;
  
/* Check for proper number of arguments */
  
   if (nrhs < 7) {
       mexErrMsgTxt("% ERR FINE: needs at least 7 input arguments.");
   } else if (nlhs < 3) {
       mexErrMsgTxt("% ERR FINE: needs 3 output argument.");
   }

/* special to set debugging requirements */

   if (nrhs <= 7) 
      idebug= 0;
   else
   {
      idebug_pt= mxGetPr(prhs[7]);
      idebug = (int)*idebug_pt;
   }
 
   if(idebug>=1) 
   {
      fprintf(stderr," nrhs = %d\n",nrhs);
      fprintf(stderr," nlhs = %d\n",nlhs);
   } 
/* ---------------------------------------------------------------------
   get vector size based on 1st input vector (assumes others same size)

      imhere= 1
*/
   m= mxGetM(prhs[0]);
   n= mxGetN(prhs[0]);
   mn= m*n;

   for(k=0;k<7; ++k)
   {
      if (mxGetM(prhs[k])!=m || mxGetN(prhs[k])!=n) {
         mexErrMsgTxt("ERROR fine1: all inputs must be same size.");
      }
   }

/*
       print*,'m,n,mn',m,n,mn
  ---------------------------------------------------------------------
     get pointer from matlab and check they are all same size
*/

   ri_pt = mxGetPr(prhs[0]);
   ro_pt = mxGetPr(prhs[1]);
   zl_pt = mxGetPr(prhs[2]);
   zu_pt = mxGetPr(prhs[3]);
   cur_pt= mxGetPr(prhs[4]);
   rp_pt = mxGetPr(prhs[5]);
   zp_pt = mxGetPr(prhs[6]);

/*
  ---------------------------------------------------------------------
   create Left Hand Side MATLAB Matrix
*/

   plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
   plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
   plhs[2] = mxCreateDoubleMatrix(m,n,mxREAL);

/*
   Get Left Hand Side Pointers
*/

   hr_pt = mxGetPr(plhs[0]);
   hz_pt = mxGetPr(plhs[1]);
   fl_pt = mxGetPr(plhs[2]);

/* 
Perform matrix vector operation in fine1.c has do loop around fine.c.
*/
   fine1(&mn,ri_pt, ro_pt, zl_pt, zu_pt, cur_pt, rp_pt, zp_pt, 
                     hr_pt, hz_pt, fl_pt);
   return;
}

/* fine1.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

int fine1(integer *num,doublereal *ri_pt, doublereal *ro_pt, doublereal *zl_pt, 
   doublereal *zu_pt, doublereal *cur_pt, doublereal *rp_pt, doublereal *zp_pt, 
   doublereal *hr_pt, doublereal *hz_pt, doublereal *fl_pt)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, m1;
    static float fl_4__, hr_4__, ri_4__, ro_4__, rp_4__, hz_4__, zl_4__, 
	    zp_4__, zu_4__;
    extern /* Subroutine */ int fine_(int *, float *, float *, float *, float 
	    *, float *, float *, float *, float *, float *, float *);
    static float cur_4__;

/*     routine to call fine many times for vectors of coils and field points */
/*     assumes the correct address is passed in arrays using %VAL() */

/*  INPUTS: */
/*      num                 = size of each vector */
/*      ri_8,ro_8,zl_8,zu_8 = inner,outer,lower,upper dimensions of coil [vector] */
/*      cur_8               = current [vector] */
/*      rp_8,zp_8           = receiving point [vector] */

/*  OUTPUTS: */
/*           hr_8, hz_8     = fields at points rp,zp [Amp/m] */
/*           fl_8           = flux at point [Vs/amu_o] */
/*                            amu_o= 4*pi*1e-7 */

/*  CALLED BY:              fine1g.f = gateway rotuine to produce MEX file */
/*  CALLS:                  fine.f = group of routines to do magnetics */

/* GA Jim Leuer 4-95 */
/* --------------------------------------------------------------------- */

/*     MATLAB real*8 (num) elements storage provided by calling routine */

/*     FORTRAN input real*4 temporary storage variables */
/*     some integer variables */
/*      integer ithru */
/* --------------------------------------------------------------------- */

/*        data ithru/0/ */
/*        ithru= ithru+1 */
/*        print*,' ' */
/*        print*,' ithru=',ithru */
/* --------------------------------------------------------------------- */
/* Main DO loop around FINE */

/*       imhere= 4 */
    /* Parameter adjustments */
    --fl_pt;
    --hz_pt;
    --hr_pt;
    --zp_pt;
    --rp_pt;
    --cur_pt;
    --zu_pt;
    --zl_pt;
    --ro_pt;
    --ri_pt;

    /* Function Body */
    i__1 = *num;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       convert from real*8 to real*4 */
	ri_4__ = ri_pt[i__];
	ro_4__ = ro_pt[i__];
	zl_4__ = zl_pt[i__];
	zu_4__ = zu_pt[i__];
	cur_4__ = cur_pt[i__];
	rp_4__ = rp_pt[i__];
	zp_4__ = zp_pt[i__];
/* 	print*,' ' */
/* 	print*,' i=',i,' ri,ro,zl,zu',ri_4,ro_4,zl_4,zu_4 */
/* 	print*,'         cur,rp,zp',cur_4,rp_4,zp_4 */
	fine_(&m1, &ri_4__, &ro_4__, &zl_4__, &zu_4__, &cur_4__, &rp_4__, &
		zp_4__, &hr_4__, &hz_4__, &fl_4__);
/* 	print*,'    m1,hr,hz,fl',m1,hr_4,hz_4,fl_4 */
/* 	print*,' ' */

/*      use existing arrays to store new data (?could dump into ri,ro,zl?) */
	hr_pt[i__] = hr_4__;
	hz_pt[i__] = hz_4__;
	fl_pt[i__] = fl_4__;
    }

    return 0;
} /* fine1_ */

