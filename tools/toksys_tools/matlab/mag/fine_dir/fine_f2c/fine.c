/*************************************************************************
FILE: fine.c

Translated from fine.f by Brian Sammuli, 04/28/08

Starting point was translation via f2c. Then manipulated by hand with the 
following changes:
	-Remove calls to f2c library functions to avoid linking to libf2c.
    -Removed extern function declarations from inside of function bodies,
    replaced with prototypes.
    -Copied d_sign function source because it is so small that it's not
    worth linking to the libf2c library.
    - Systematically removed all goto's and tested with test_finec.m.
    - Reindexed for loops.

****************************************************************************/

#include "math.h"
#include "f2c.h"

/*----------------------------------------------------------------------------
Common Block Declarations 
----------------------------------------------------------------------------*/

struct xw_1_ {
    real x[182], w[182];
    integer n1[40];
};

struct xw_2_ {
    real x1[86], x2[96], w1[86], w2[96];
    integer n1[40];
};

struct qep_1_ {
    doublereal aq[12], ae[12], ap[12];
};


/*----------------------------------------------------------------------------
Defines
----------------------------------------------------------------------------*/
#define DEBUGGING 0

#define xw_1 (*(struct xw_1_ *) &xw_)
#define xw_2 (*(struct xw_2_ *) &xw_)
#define qep_1 (*(struct qep_1_ *) &qep_)


/*----------------------------------------------------------------------------
Function Prototypes
----------------------------------------------------------------------------*/

double d_sign(doublereal *a, doublereal *b);
int ohfiga_(integer *, real *, real *);
int ohfilf_(real *, real *, real *, real *, real *, real *, real *);


/*----------------------------------------------------------------------------
Initialized data 
----------------------------------------------------------------------------*/
struct {
    real e_1[364];
    integer e_2[40];
    } xw_ = { .5773503f, -.5773503f, .7745967f, 0.f, -.7745967f, .8611363f, 
	    .339981f, -.339981f, -.8611363f, .9061798f, .5384693f, 0.f, 
	    -.5384693f, -.9061798f, .9324695f, .6612094f, .2386192f, 
	    -.2386192f, -.6612094f, -.9324695f, .9602899f, .7966665f, 
	    .5255324f, .1834346f, -.1834346f, -.5255324f, -.7966665f, 
	    -.9602899f, .9739065f, .8650634f, .6794096f, .4333954f, .1488743f,
	     -.1488743f, -.4333954f, -.6794096f, -.8650634f, -.9739065f, 
	    .9815606f, .9041173f, .7699027f, .587318f, .3678315f, .1252334f, 
	    -.1252334f, -.3678315f, -.587318f, -.7699027f, -.9041173f, 
	    -.9815606f, .9894009f, .944575f, .8656312f, .7554044f, .6178762f, 
	    .4580168f, .2816036f, .0950125f, -.0950125f, -.2816036f, 
	    -.4580168f, -.6178762f, -.7554044f, -.8656312f, -.944575f, 
	    -.9894009f, .9931286f, .9639719f, .9122344f, .839117f, .7463319f, 
	    .6360537f, .510867f, .3737061f, .2277859f, .0765265f, -.0765265f, 
	    -.2277859f, -.3737061f, -.510867f, -.6360537f, -.7463319f, 
	    -.839117f, -.9122344f, -.9639719f, -.9931286f, .9951872f, 
	    .9747286f, .9382746f, .8864155f, .820002f, .7401242f, .6480937f, 
	    .5454215f, .4337935f, .3150427f, .1911189f, .0640569f, -.0640569f,
	     -.1911189f, -.3150427f, -.4337935f, -.5454215f, -.6480937f, 
	    -.7401242f, -.820002f, -.8864155f, -.9382746f, -.9747286f, 
	    -.9951872f, .9972639f, .9856115f, .9647623f, .9349061f, .8963212f,
	     .8493676f, .7944838f, .7321821f, .6630443f, .5877158f, .5068999f,
	     .4213513f, .3318686f, .2392874f, .144472f, .0483077f, -.0483077f,
	     -.144472f, -.2392874f, -.3318686f, -.4213513f, -.5068999f, 
	    -.5877158f, -.6630443f, -.7321821f, -.7944838f, -.8493676f, 
	    -.8963212f, -.9349061f, -.9647623f, -.9856115f, -.9972639f, 
	    .9982377f, .9907262f, .9772599f, .9579168f, .9328128f, .9020988f, 
	    .8659595f, .8246122f, .7783057f, .7273183f, .6719567f, .6125539f, 
	    .5494671f, .4830758f, .4137792f, .3419941f, .2681522f, .1926976f, 
	    .1160841f, .0387724f, -.0387724f, -.1160841f, -.1926976f, 
	    -.2681522f, -.3419941f, -.4137792f, -.4830758f, -.5494671f, 
	    -.6125539f, -.6719567f, -.7273183f, -.7783057f, -.8246122f, 
	    -.8659595f, -.9020988f, -.9328128f, -.9579168f, -.9772599f, 
	    -.9907262f, -.9982377f, 1.f, 1.f, .5555556f, .8888889f, .5555556f,
	     .3478548f, .6521452f, .6521452f, .3478548f, .2369269f, .4786287f,
	     .5688889f, .4786287f, .2369269f, .1713245f, .3607616f, .4679139f,
	     .4679139f, .3607616f, .1713245f, .1012285f, .222381f, .3137066f, 
	    .3626838f, .3626838f, .3137066f, .222381f, .1012285f, .0666713f, 
	    .1494513f, .2190864f, .2692667f, .2955242f, .2955242f, .2692667f, 
	    .2190864f, .1494513f, .0666713f, .0471753f, .1069393f, .1600783f, 
	    .2031674f, .2334925f, .249147f, .249147f, .2334925f, .2031674f, 
	    .1600783f, .1069393f, .0471753f, .0271525f, .0622535f, .0951585f, 
	    .124629f, .149596f, .1691565f, .1826034f, .1894506f, .1894506f, 
	    .1826034f, .1691565f, .149596f, .124629f, .0951585f, .0622535f, 
	    .0271525f, .017614f, .0406014f, .062672f, .0832767f, .1019301f, 
	    .1181945f, .1316886f, .1420961f, .149173f, .1527534f, .1527534f, 
	    .149173f, .1420961f, .1316886f, .1181945f, .1019301f, .0832767f, 
	    .062672f, .0406014f, .017614f, .0123412f, .0285314f, .0442774f, 
	    .0592986f, .0733465f, .0861902f, .0976187f, .1074443f, .1155057f, 
	    .1216705f, .1258375f, .1279382f, .1279382f, .1258375f, .1216705f, 
	    .1155057f, .1074443f, .0976187f, .0861902f, .0733465f, .0592986f, 
	    .0442774f, .0285314f, .0123412f, .0070186f, .0162744f, .0253921f, 
	    .0342739f, .0428359f, .0509981f, .0586841f, .0658222f, .0723458f, 
	    .0781939f, .0833119f, .0876521f, .0911739f, .0938444f, .0956387f, 
	    .0965401f, .0965401f, .0956387f, .0938444f, .0911739f, .0876521f, 
	    .0833119f, .0781939f, .0723458f, .0658222f, .0586841f, .0509981f, 
	    .0428359f, .0342739f, .0253921f, .0162744f, .0070186f, .0045213f, 
	    .0104983f, .0164211f, .0222458f, .027937f, .0334602f, .0387822f, 
	    .0438709f, .0486958f, .0532278f, .0574398f, .0613062f, .064804f, 
	    .067912f, .0706116f, .0728866f, .0747232f, .0761104f, .0770398f, 
	    .0775059f, .0775059f, .0770398f, .0761104f, .0747232f, .0728866f, 
	    .0706116f, .067912f, .064804f, .0613062f, .0574398f, .0532278f, 
	    .0486958f, .0438709f, .0387822f, .0334602f, .027937f, .0222458f, 
	    .0164211f, .0104983f, .0045213f, 500, 0, 2, 5, 9, 14, 500, 20, 
	    500, 28, 500, 38, 500, 500, 500, 26, 500, 500, 500, 66, 500, 500, 
	    500, 86, 500, 500, 500, 500, 500, 500, 500, 110, 500, 500, 500, 
	    500, 500, 500, 500, 142 };

struct {
    doublereal e_1[36];
    } qep_ = { 3., -.375, -.046875, -.0146484375, -.00640869140625, 
	    -.00336456298828125, -.00198268890380859, -.00126573443412781, 
	    -8.57007689774036e-4, -6.07047113589943e-4, -4.45627767476255e-4, 
	    -3.36752801558762e-4, .25, .015625, .00390625, .00152587890625, 
	    7.476806640625e-4, 4.20570373535156e-4, 2.59637832641601e-4, 
	    1.71401537954807e-4, 1.19028845801949e-4, 8.59983410919084e-5, 
	    6.4143390773097e-5, 4.91097835606524e-5, 1., .375, .234375, 
	    .1708984375, .13458251953125, .111030578613281, .0945081710815427,
	     .0822727382183072, .0728456536307928, .0653587392298502, 
	    .0592684930743416, .0542172010509602 };


/*----------------------------------------------------------------------------
Table of constant values 
----------------------------------------------------------------------------*/

static real c_b16 = -1.1e4f;
static real c_b17 = 1.1e4f;
static integer c__1 = 1;
static doublereal c_b48 = .25;
static doublereal c_b49 = .5;



/*----------------------------------------------------------------------------
Function Definitions
----------------------------------------------------------------------------*/

/*
Function: fine_
*/

int fine_(integer *m1, real *a3, real *a4, real *z3, real *
	z4, real *at1, real *rp, real *zp, real *hr, real *hz, real *fl)
{
#if 0
    /* Format strings */
    static char fmt_31[] = "(\002 OR IS LESS THAN IR (A1 = \002,f6.3,\002A2 "
	    "= \002,f6.3,\002)\002)";
    static char fmt_33[] = "(\002 EXIT IN SUBR. OHFINE \002,6e12.5)";
#endif

    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static real a;
    static integer j, n;
    static real q, a1, a2, z1, z2, ac, da, cd, at, pi;
    static integer it;
    static real rq, sq, st, zq, fl1[40], hr1[40], hz1[40], fll, fln, flp, hrl,
	     hrn, hrp, hzl, hzn, hzp;
    static logical inside;

    #if 0 
    /* Fortran I/O blocks */
    static cilist io___32 = { 0, 6, 0, fmt_31, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_33, 0 };
    #endif
    
    ohfidt_();
    pi = 3.141593f;
    j = 3;
    a1 = *a3;
    a2 = *a4;
    at = *at1;
    inside = FALSE_;
    if (!(*a3 >= *rp || *rp >= *a4)) {
    /*     IS POINT INSIDE COIL ? */
        if (*z3 < *zp && *zp < *z4) {
            inside = TRUE_;
        }
    /*     SPLIT COIL AT RP */
    /*     CAL. LOWER PORTION */
        j = 1;
        at = *at1 * (*rp - *a3) / (*a4 - *a3);
        a1 = *a3;
        a2 = *rp;
    }
/*     STM. 3 TO 20 CONSTITUTES A SINGLE COIL CASE */

    while(1) {
        if (*z4 < 1e4f) {
            if (!(*z3 > -1e4f)) {
                z1 = -1e4f;
                z2 = *z4 - *zp;
                zq = dabs(z2);
                cd = at;             
            } else {
                z1 = *z3 - *zp;
                z2 = *z4 - *zp;
                if (*z3 != *z4) {
                    cd = at / (*z4 - *z3);
                }
                
                /*     DETERMINE GAUSS QUART. ORDER M1 */
                zq = dabs(z1);
                if (*zp - (*z3 + *z4) * .5f > 0.f) {
                    zq = dabs(z2);
                }
            }
            
        } else {
            z1 = *z3 - *zp;
            z2 = 1e4f;
            zq = dabs(z1);
            cd = at;
        }
        ac = (a2 + a1) / 2.f;
        da = a2 - a1 + 1e-10f;
        if (da < 0.f) {
            return 0;
        }
        rq = (r__1 = *rp - ac, dabs(r__1)) - da / 2.f;
        q = dmax(rq,zq);
        if (q < 1e-5f) {
            q = 1e-5f;
        }
        sq = da / q;
        if (sq > 4.f) {
            sq = 8.f;
        }
        st = sq * da / ac;
        *m1 = 3;
        if ((0.f < st && st <= .1f)) {
            *m1 = 3;
        } else if ((.1f < st && st <= .4f)) {
            *m1 = 4;
        } else if ((.4f < st && st <= 1.f)) {
            *m1 = 5;
        } else if ((1.f < st && st <= 2.f)) {
            *m1 = 6;
        } else if ((2.f < st && st <= 4.f)) {
            *m1 = 8;
        } else if ((4.f < st && st <= 8.f)) {
            *m1 = 10;
        } else if ((8.f < st && st <= 1e9f)) {
            *m1 = 12;   
        } else {
            return 0;
        }

        if (!inside) {
        /*     POINT OUTSIDE SOURCE REGION. */
            i__1 = *m1;
            for (n = 0; n < i__1; ++n) {
                a = da * .5f * (xw_1.x[xw_1.n1[*m1 - 1] + n] + 1.f) + a1;
                if (z1 != z2) {
                    ohfisl_(&a, &z1, &z2, &cd, rp, &hr1[n], &hz1[n], &fl1[n], 
                            &it);
                } else {
                    r__1 = -z1;
                    ohfilf_(&a, &r__1, &at, rp, &hr1[n], &hz1[n], &fl1[n]);
                }
            }
        } else {
        /*     POINT INSIDE SOURCE REGION. */
            i__1 = *m1;
            for (n = 0; n < i__1; ++n) {
                a = da * .5f * (xw_1.x[xw_1.n1[*m1 - 1] + n] + 1.f) + a1;
                ohfisl_(&a, &c_b16, &z1, &cd, rp, &hrn, &hzn, &fln, &it);
                ohfisl_(&a, &z2, &c_b17, &cd, rp, &hrp, &hzp, &flp, &it);
                hr1[n] = hrn + hrp;
                hz1[n] = hzn + hzp;
                fl1[n] = fln + flp;
            }
        }
        ohfiga_(m1, hr1, hr);
        ohfiga_(m1, hz1, hz);
        ohfiga_(m1, fl1, fl);

        if (j == 1) {
            /*     CAL. UPPER PORTION */
            hrl = *hr;
            hzl = *hz;
            fll = *fl;
            j = 2;
            at = *at1 * (*a4 - *rp) / (*a4 - *a3);
            a1 = *rp;
            a2 = *a4;
        } else if (j == 2) {
            if (inside) {
                cd = *at1 / ((*z4 - *z3) * (*a4 - *a3));
                *hr = -(*hr) - hrl;
                *hz = cd * (*a4 - *rp) - *hz - hzl;
                /* Computing 2nd power */
                r__1 = *rp;
                /* Computing 3rd power */
                r__2 = *a3;
                *fl = cd * pi * (r__1 * r__1 * (*a4 * 3.f - *rp * 2.f) - r__2 * (r__2 * 
                    r__2)) / 3.f - *fl - fll;
                
            } else { /* outside */
                *hr += hrl;
                *hz += hzl;
                *fl += fll;
            } 
            return 0;    
        } else if (j == 3) {
            return 0;
        } else { /* j has assumed an illegal value */
            return 1;
        }
    }
} /* fine_ */




/*
Function: ohfidt_
*/
int ohfidt_(void)
{
/*     BLOCK DATA */
/*     GAUSSIAN QUAD. ORDERS, - 2,3,4,5,6,8,10,12,16,20,24,32,40 */
/*     COEFFICIENTS FOR FILAMENT FIELD. */
    return 0;
} /* ohfidt_ */



/*
Function: ohfiga_
*/
int ohfiga_(integer *m1, real *y, real *int__)
{
#if 0
    /* Format strings */
    static char fmt_1[] = "(\002  M1 = \002,i5,\002NOT STANDARD\002)";
#endif

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

    #if 0
    /* Fortran I/O blocks */
    static cilist io___35 = { 0, 6, 0, fmt_1, 0 };
    #endif

    /* Parameter adjustments */
    --y;

    /* Function Body */

    if (xw_1.n1[*m1 - 1] == 500) {
	    return 0;
    }

    *int__ = 0.f;
    i__1 = *m1;
    for (n = 0; n < i__1; ++n) {

	    *int__ += xw_1.w[xw_1.n1[*m1 - 1] + n] * y[n + 1];
    }
    *int__ *= .5f;
    return 0;
} /* ohfiga_ */



/*
Function: ohfisl_
*/
int ohfisl_(real *as, real *z1, real *z2, real *cds, real *
	rs, real *hr, real *hz, real *fl, integer *it)
{
#if 0
    /* Format strings */
    static char fmt_1[] = "(//\002   POINT AT EDGE IN OHFISL, A,Z1,Z2,Z,R,R1"
	    "2,APR2,K2 \002,/(2x,6d21.15))";
#endif

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal a;
    static integer m;
    static doublereal r__, s, z__, k2, r1, cd, r12, pi, fl1, fl2, hr1, hr2, 
	    hz1, hz2, fac, alf, del, gam, bet, amr, eps, zet, alf1, amr2, 
	    apr2, zet1;

    #if 0
    /* Fortran I/O blocks */
    static cilist io___64 = { 0, 6, 0, fmt_1, 0 };
    #endif

    /*printf("Inside ohfisl_\n");*/

/*     9/6/73 REVISED VERSION   U. R. CHRISTENSEN */
    a = *as;
    cd = *cds;
    r__ = *rs;
    pi = 3.141592653589793;
    m = 1;
    z__ = *z2;

    amr = a - r__;
    if (abs(amr) < 5e-7) {
	amr = 0.;
    }
    /* Computing 2nd power */
    d__1 = amr;
    amr2 = d__1 * d__1;
    /* Computing 2nd power */
    d__1 = a + r__;
    apr2 = d__1 * d__1;


/*     THE FOLLOWING STATEMENTS ARE TO SPEED-UP THE SEMI-INFINITE CASE. */

    while(1) {
        if (abs(z__) < 1e4f) {
            r12 = apr2 + z__ * z__;
            r1 = sqrt(r12);
            k2 = a * 4. * r__ / r12;
            hr1 = 0.;
            hz1 = 0.;
            fl1 = 0.;
            *it = 0;
            if (!(k2 >= 1.)) {
                /* POINT AT EDGE,  NO CONTRIBUTION. */

                alf = 1.;
                bet = sqrt(1. - k2);
                del = 1.;
                if (amr2 != 0.) {
                    del = amr2 / apr2 / bet;
                    eps = a * 4. * r__ / amr2;
                }
                zet = 0.;
                s = 0.;
                fac = .5;
                
                while(1) {
                    /* Computing 2nd power */
                    d__1 = alf - bet;
                    s += fac * (d__1 * d__1);
                    if (*it >= 20) {
                        break;
                    }
                    alf1 = alf;
                    alf = (alf + bet) * .5;
                    bet = sqrt(alf1 * bet);
                    zet1 = zet;
                    zet = (eps + zet) * .5;
                    eps = (del * eps + zet1) / (del + 1.);
                    del = (del + 2. + 1. / del) * bet / (alf * 4.);
                    fac *= 2.;
                    ++(*it);
                    if ((d__1 = alf - bet, abs(d__1)) > 1e-7f) {
                        continue;
                    }
                    if (z__ == 0.) {
                        break;
                    }
                    if ((d__1 = del - 1., abs(d__1)) > 1e-7f) {
                        continue;
                    }
                    break;
                }
                gam = s + k2;
                hr1 = 0.;
                if (r__ != 0.f) {
                    hr1 = r1 * s / (r__ * 8. * alf);
                }
                hz1 = z__ * (a * 2. + amr * zet) / ((a + r__) * 4. * r1 * alf);
                fl1 = pi * z__ * (r1 * .5 * gam - amr2 * zet / r1) / (alf * 4.);
            
            }    

        } else {
            hr1 = 0.;
        
            if (amr < 0.) {
                hz1 = 0.;
                d__1 = pi * .5 * a * a;
                fl1 = d_sign(&d__1, &z__);

            } else if (amr == 0) {
                hz1 = d_sign(&c_b48, &z__);
                d__1 = pi * .5 * a * a;
                fl1 = d_sign(&d__1, &z__);

            } else {
                hz1 = d_sign(&c_b49, &z__);
                d__1 = pi * .5 * r__ * r__;
                fl1 = d_sign(&d__1, &z__);

            }
        }
        
        if (m != 1) {
            *hr = cd * (hr2 - hr1);
            *hz = cd * (hz2 - hz1);
            *fl = cd * (fl2 - fl1);
            return 0;
        }
        m = 2;
        z__ = *z1;
        hr2 = hr1;
        hz2 = hz1;
        fl2 = fl1;
    }
} /* ohfisl_ */




/* 
Function: ohfilf_
 */ 
int ohfilf_(real *a, real *z__, real *ni, real *r__, real *
	hr, real *hz, real *fl)
{
    /* Initialized data */

    static doublereal ak1 = 1.38629436112;
    static doublereal ak2 = .09666344259;
    static doublereal ak3 = .03590092383;
    static doublereal ak4 = .03742563713;
    static doublereal ak5 = .01451196212;
    static doublereal bk1 = .5;
    static doublereal bk2 = .12498593597;
    static doublereal bk3 = .06880248576;
    static doublereal bk4 = .03328355346;
    static doublereal bk5 = .00441787012;
    static doublereal ae1 = .44325141463;
    static doublereal ae2 = .0626060122;
    static doublereal ae3 = .04757383546;
    static doublereal ae4 = .01736506451;
    static doublereal be1 = .2499836831;
    static doublereal be2 = .09200180037;
    static doublereal be3 = .04069697526;
    static doublereal be4 = .00526449639;

    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real b;
    static doublereal e, k;
    static real m;
    static integer n;
    static real p, q;
    static doublereal v, w;
    static real c2, m2, r1;
    static doublereal ea, eb, ka, kb;
    static real ke, r12, r22, re, pu, qt, err;

    err = 1e-6f;
    if (*ni == 0.f) {
	*hr = 0.f;
        *hz = 0.f;
        *fl = 0.f;
        return 0;
    }
/* Computing 2nd power */
    r__1 = *a + *r__;
/* Computing 2nd power */
    r__2 = *z__;
    r12 = r__1 * r__1 + r__2 * r__2;
    r1 = sqrt(r12);
    r22 = r12 - *a * 4.f * *r__;
    c2 = *a * 4.f * *r__ / r12;
    if (c2 >= 1.f) {  
        *hr = 1e10f;
        *hz = 1e10f;
        *fl = 1e10f;
        return 0;
    }

    if (c2 < .9f) {
        n = 0;
        /* Computing 2nd power */
        r__1 = r1 + sqrt((dabs(r22)));
        m = *a * 4.f * *r__ / (r__1 * r__1);
        m2 = m * m;
        pu = 1.f;
        q = 0.f;
        e = 1.f;
        p = 0.f;
        b = *ni / (r1 * 4.f * (1.f - m) * (1.f - m2));
        *hr = 0.f;
        if (*r__ != 0.f) {
            for (n = 0; n < 12; ++n) {
                pu *= m2;
                qt = qep_1.aq[n] * pu;
                q += qt;
                e += qep_1.ae[n] * pu;
                p += qep_1.ap[n] * pu;
                if ((r__1 = qt / q, dabs(r__1)) < err) {
                    break;
                }
            }
	        *hr = b * *z__ * q / *r__;
	    }
	/* Computing 2nd power */
        r__1 = *a * (m + 1.f);
        *hz = b * (e * 2.f * (r__1 * r__1) / r12 - q);
        *fl = *ni * 1.5707963f * r1 * p / (m + 1.f);
        return 0;

    }
    v = (doublereal) (r22 / r12);
    ea = (((ae4 * v + ae3) * v + ae2) * v + ae1) * v + 1.;
    eb = (((be4 * v + be3) * v + be2) * v + be1) * v;
    ka = (((ak5 * v + ak4) * v + ak3) * v + ak2) * v + ak1;
    kb = (((bk5 * v + bk4) * v + bk3) * v + bk2) * v + bk1;
    w = log(v);
    e = ea - w * eb;
    k = ka - w * kb;
    ke = (real) (k - e);
    re = *a * 2.f * (real) e / r22;
    b = *ni * .15915494f / r1;
    *hr = 0.f;
    if (*r__ != 0.f) {
	    *hr = b * (-ke + *r__ * re) * *z__ / *r__;
    }

    *hz = b * (ke + (*a - *r__) * re);
    *fl = *ni * r1 * (ke - c2 * .5f * (real) k);
    return 0;

} /* ohfilf_ */


/*
Function: d_sign
*/
double d_sign(doublereal *a, doublereal *b)
{
    double x;
    
    x = (*a >= 0 ? *a : - *a);
    return( *b >= 0 ? x : -x);
}

