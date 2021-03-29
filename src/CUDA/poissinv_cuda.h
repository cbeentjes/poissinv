/*

    This software was written by Mike Giles, copyright University of Oxford,
    and is provided under the terms of the GNU GPLv3 license:
    http://www.gnu.org/licenses/gpl.html

    Commercial users who would like to use the software under a more
    permissive license, such as BSD, should contact the author:
    mike.giles@maths.ox.ac.uk

*/

// include CUDA header which defines Inf and NaN constants

#include <math_constants.h>


//
// This single precision function computes the inverse
// of the Poisson CDF, using at most 10 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~1.6e+07
//
// Compilation with use_fast_math. Extra error for 5.5 < lam <= 10 due to
// conversion of epxf(0.5f*lam) into single precision intrinsics
// __expf(0.5f*lam)in the direct summation.
// ave |error| < 6.2e-07*sqrt(max(1e-02,lam)   for 0   <= lam <= 5.5
//             < 2.0e-06                       for 5.5 <  lam <= 10
//             < 1.8e-07*sqrt(max(69   ,lam)   for 10  <  lam <  ~1.6e+07
//
// Compilation without use_fast_math
// ave |error| < 6.2e-07*sqrt(max(1e-02,lam)   for 0    <= lam <= 7.5
//             < 1.8e-07*sqrt(max(89   ,lam)   for 7.5  <  lam <  ~1.6e+07
//
// For lam > ~1.6e+07, the errors will be about 1 ulp.
//

__device__ static double rlog1(double);

__device__ inline float poissinvf(float u, float lam) {

    float s, t, x=0.0f;

    // handle exceptions

    if (lam < 0.0f || isinf(lam)) return CUDART_NAN_F;
    if (u  < 0.0f || u > 1.0f) return CUDART_NAN_F;
    if (u == 0.0f || lam == 0.0f) return 0.0f;
    if (u == 1.0f) return CUDART_INF_F;

    // large lam
    // NOTE: need threshold for large Lam > 4.0.
    // For fixed lam the threshold Lam > 10.0 is faster
    // but for mixed lam it is better to choose Lam > 4.0 due
    // to warp divergence

    if (lam > 4.0f) {
        s = normcdfinvf(u)*rsqrtf(lam);

        // use polynomial approximations in central region

        if ((s > -0.6833501f) && (s < 1.777993f)) {
            float rm;

            //  polynomial approximation to f^{-1}(s) - 1

            rm =  2.82298751e-07f;
            rm = -2.58136133e-06f + rm*s;
            rm =  1.02118025e-05f + rm*s;
            rm = -2.37996199e-05f + rm*s;
            rm =  4.05347462e-05f + rm*s;
            rm = -6.63730967e-05f + rm*s;
            rm =  0.000124762566f + rm*s;
            rm = -0.000256970731f + rm*s;
            rm =  0.000558953132f + rm*s;
            rm =  -0.00133129194f + rm*s;
            rm =   0.00370367937f + rm*s;
            rm =   -0.0138888706f + rm*s;
            rm =     0.166666667f + rm*s;
            rm =             s + s*(rm*s);

            //  polynomial approximation to correction c0(r)

            t  =   1.86386867e-05f;
            t  =  -0.000207319499f + t*rm;
            t  =     0.0009689451f + t*rm;
            t  =   -0.00247340054f + t*rm;
            t  =    0.00379952985f + t*rm;
            t  =   -0.00386717047f + t*rm;
            t  =    0.00346960934f + t*rm;
            t  =   -0.00414125511f + t*rm;
            t  =    0.00586752093f + t*rm;
            t  =   -0.00838583787f + t*rm;
            t  =     0.0132793933f + t*rm;
            t  =     -0.027775536f + t*rm;
            t  =      0.333333333f + t*rm;

            //  O(1/lam) correction

            x  =  -0.000134529865f;
            x  =     0.0013711034f + x*rm;
            x  =   -0.00583140335f + x*rm;
            x  =      0.013452415f + x*rm;
            x  =    -0.0185780057f + x*rm;
            x  =     0.0169826221f + x*rm;
            x  =    -0.0135307848f + x*rm;
            x  =     0.0135629517f + x*rm;
            x  =    -0.0155146497f + x*rm;
            x  =     0.0174051522f + x*rm;
            x  =    -0.0198016963f + x*rm;
            x  = __fdividef(x,lam);

            //    sum from smallest to largest to minimise rounding error;
            //    use of __fadd_rd to round down final sum is important for
            //    very large values of lambda to ensure correct rounding

            x = truncf( __fadd_rd(lam, (x+t)+lam*rm) );
        }

        // otherwise use Newton iteration

        else if (s > -sqrtf(2.0f)) {
            float r, r2, s2;

            r = 1.0f + s;
            if (r < 0.1f) r = 0.1f;

            do {
                t  = __logf(r);
                r2 = r;
                s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                if (r < 1.0f) s2 = -s2;
                r = r2 - (s2-s)*s2/t;
                if (r < 0.1f*r2) r = 0.1f*r2;
            } while (fabsf(r-r2) > 1e-5f);

            t = __logf(r);
            x = lam*r + __logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;
            x = truncf( x - 0.0218f/(x+0.065f*lam) );
        }
    }

    // bottom-up summation

    if (x < 10.0f) {
        float del;

        x    = 0.0f;
        t    = expf(0.5f*lam);
        del  = 0.0f;
        if (u > 0.5f)
            del  = t*(1e-6f*t);
        s    = 1.0f - t*(u*t) + del;

        while (s < 0.0f) {
            x  += 1.0f;
            t   = x/lam;
            del = t*del;
            s   = t*s + 1.0f;
        }

        // top-down summation if needed

        if (s < 2.0f*del) {
            del = 1e6f*del;
            t   = 1e7f*del;
            del = (1.0f-u)*del;

            while (del < t) {
                x   += 1.0f;
                del *= x/lam;
            }

            s = del;
            t = 1.0f;
            while (s > 0.0f) {
                t *= x/lam;
                s -= t;
                x -= 1.0f;
            }
        }
    }

    return x;
}


//
// This double precision function computes the inverse
// of the Poisson CDF, using about 30 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~9e+15
// ave |error| < 1.2e-15*sqrt(max(1e-02,lam))    for 0 <= lam <= 6
//             < 2.0e-16*sqrt(max(216  ,lam))    for 6 <  lam < ~9e+15
//
// For lam >= ~9e+15, the errors will be about 1 ulp.
//

//
// naming convention: double precision variables are capitalised
//

__device__ inline double poissinv(double U, double Lam) {

    double X=0.0, Xi, S, T, Rm;

    // handle exceptions

    if ((Lam < 0.0) || isinf(Lam)) return CUDART_NAN;
    if (U <  0.0 || U > 1.0) return CUDART_NAN;
    if (U == 0.0 || Lam == 0.0) return 0.0;
    if (U == 1.0) return CUDART_INF;

    // large lam
    // NOTE: need threshold for large Lam > 4.0.
    // For fixed lam the threshold Lam > 15.0 is faster
    // but for mixed lam it is better to choose Lam > 4.0 due
    // to warp divergence

    if (Lam > 4.0) {
        float s, t, del, rm;

        S   = normcdfinv(U)*rsqrt(Lam);
        s   = (float) S;

        // del = 2.095e-06f is sufficient for U > 1e-195, but for very small U
        // a larger del = 3.3e-06f is correct. Minimal effect on samples/sec
        // from changing del to the smaller of the two.

        del = 3.3e-6f;

        // use polynomial approximations in central region

        if ((s > -0.6833501f) && (s < 1.777993f)) {
            float x;

            //  polynomial approximation to f^{-1}(s) - 1

            rm =    2.82298751e-07f;
            rm =   -2.58136133e-06f + rm*s;
            rm =    1.02118025e-05f + rm*s;
            rm =   -2.37996199e-05f + rm*s;
            rm =    4.05347462e-05f + rm*s;
            rm =   -6.63730967e-05f + rm*s;
            rm =    0.000124762566f + rm*s;
            rm =   -0.000256970731f + rm*s;
            rm =    0.000558953132f + rm*s;
            rm =    -0.00133129194f + rm*s;
            rm =     0.00370367937f + rm*s;
            rm =     -0.0138888706f + rm*s;
            Rm = 0.1666666666666667 + rm*s;

            S +=                  S*(Rm*S);
            Rm = S;
            rm = (float) Rm;

            //  polynomial approximation to correction c0(r)

            t  =   1.86386867e-05f;
            t  =  -0.000207319499f + t*rm;
            t  =     0.0009689451f + t*rm;
            t  =   -0.00247340054f + t*rm;
            t  =    0.00379952985f + t*rm;
            t  =   -0.00386717047f + t*rm;
            t  =    0.00346960934f + t*rm;
            t  =   -0.00414125511f + t*rm;
            t  =    0.00586752093f + t*rm;
            t  =   -0.00838583787f + t*rm;
            t  =     0.0132793933f + t*rm;
            t  =     -0.027775536f + t*rm;
            t  =      0.333333333f + t*rm;

            //  O(1/lam) correction

            x  =  -0.000134549989f;
            x  =    0.00137126732f + x*rm;
            x  =    -0.0058318548f + x*rm;
            x  =     0.0134526835f + x*rm;
            x  =    -0.0185770659f + x*rm;
            x  =     0.0169808982f + x*rm;
            x  =    -0.0135302102f + x*rm;
            x  =     0.0135636388f + x*rm;
            x  =      -0.01551508f + x*rm;
            x  =     0.0174051003f + x*rm;
            x  =    -0.0198016502f + x*rm;

            x  = __fdividef(x,(float) Lam);

            //    sum from smallest to largest to minimise rounding error;
            S = ((x + del) + t) + Lam*Rm;
        }

        // otherwise use Newton iteration

        else if (s > -sqrtf(2.0f)) {
            float r, r2, s2;

            r = 1.0f + s;
            if (r < 0.1f) r = 0.1f;

            do {
                t  = __logf(r);
                r2 = r;
                s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                if (r < 1.0f) s2 = -s2;
                r = r2 - (s2-s)*s2/t;
                if (r < 0.1f*r2) r = 0.1f*r2;
            } while (fabsf(r-r2) > 1e-5f);

            t   = __logf(r);
            rm  = r - 1.0f;
            t   = __logf(sqrtf(2.0f*r*(-rm+r*t))/fabsf(rm)) / t;

            // O(1/lam) ad-hoc correction
            t  += -0.0218/(Lam*(0.065 + r));
            del =  0.01/(Lam*r);

            S = (del + t) + Lam*rm;
        }

        // if x>10, round down to nearest integer, and check accuracy

        // NOTE: use of __dadd_rd to round down final sum is important
        // for large values of lambda to ensure fast execution and accuracy
        // If Lam >> 1 then del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        X = __dadd_rd(Lam, S);
        X = trunc(X);

        if (X > 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0f*del >= X - Lam) return X;

            // correction procedure based on Temme approximation (double precision)

            if (X > 0.5*Lam && X < 2.0*Lam) {
                double Eta, B0, B1;

                Xi = 1.0 / X;
                // Use rlog1(z) = z - log(1+z) function, a fast rational
                // minimax approximation for -0.6 < z < 1.6. This is a lot
                // more accurate when x is near 0 (here X near Lambda).
                Eta = sqrt(2.0*rlog1((Lam-X)*Xi));
                if (X > Lam) Eta = -Eta;

                B1 =  8.0995211567045583e-16;              S = B1;
                B0 = -1.9752288294349411e-15;              S = B0 + S*Eta;
                B1 = -5.1391118342426808e-16 + 25.0*B1*Xi; S = B1 + S*Eta;
                B0 =  2.8534893807047458e-14 + 24.0*B0*Xi; S = B0 + S*Eta;
                B1 = -1.3923887224181616e-13 + 23.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.3717632624009806e-13 + 22.0*B0*Xi; S = B0 + S*Eta;
                B1 =  1.1004392031956284e-13 + 21.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.0276692801141763e-12 + 20.0*B0*Xi; S = B0 + S*Eta;
                B1 =  2.4361948020667402e-11 + 19.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.8307721325504166e-11 + 18.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.5514193994946487e-11 + 17.0*B1*Xi; S = B1 + S*Eta;
                B0 =  9.1476995822367933e-10 + 16.0*B0*Xi; S = B0 + S*Eta;
                B1 = -4.3820360184533521e-09 + 15.0*B1*Xi; S = B1 + S*Eta;
                B0 =  1.0261809784240299e-08 + 14.0*B0*Xi; S = B0 + S*Eta;
                B1 =  6.7078535434015332e-09 + 13.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.7665952736826086e-07 + 12.0*B0*Xi; S = B0 + S*Eta;
                B1 =  8.2967113409530833e-07 + 11.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.8540622107151585e-06 + 10.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.1854485106799979e-06 +  9.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.9192631785224383e-05 +  8.0*B0*Xi; S = B0 + S*Eta;
                B1 = -0.00017875514403292177 +  7.0*B1*Xi; S = B1 + S*Eta;
                B0 =  0.00035273368606701921 +  6.0*B0*Xi; S = B0 + S*Eta;
                B1 =   0.0011574074074074078 +  5.0*B1*Xi; S = B1 + S*Eta;
                B0 =   -0.014814814814814815 +  4.0*B0*Xi; S = B0 + S*Eta;
                B1 =    0.083333333333333329 +  3.0*B1*Xi; S = B1 + S*Eta;
                B0 =    -0.33333333333333331 +  2.0*B0*Xi; S = B0 + S*Eta;
                S  = S / (1.0 + B1*Xi);

                S = S*exp(-0.5*X*Eta*Eta)*rsqrt(2.0*3.141592653589793*X);
                if (X < Lam) {
                    S += 0.5*erfc(Eta*sqrt(0.5*X));
                    if (S > U) X -= 1.0;
                }
                else {
                    S -= 0.5*erfc(-Eta*sqrt(0.5*X));
                    if (S > U-1.0) X -= 1.0;
                }
            }

            // sum downwards or upwards

            else {
                Xi = 1.0 / X;
                S = - (691.0/360360.0);
                S =   (1.0/1188.0) + S*Xi*Xi;
                S = - (1.0/1680.0) + S*Xi*Xi;
                S =   (1.0/1260.0) + S*Xi*Xi;
                S = - (1.0/360.0)  + S*Xi*Xi;
                S =   (1.0/12.0)   + S*Xi*Xi;
                S =                  S*Xi;
                S = (X - Lam) - X*log(X/Lam) - S;

                if (X < Lam) {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
                    T  = 1.0;
                    Xi = X;
                    for (int i=1; i<50; i++) {
                        Xi -= 1.0;
                        T  *= Xi/Lam;
                        S  += T;
                    }
                    if (S > 0.0) X -= 1.0;
                }

                else {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*((1.0-U)*T) * sqrt(2.0*3.141592653589793*X);
                    Xi = X;
                    for (int i=0; i<50; i++) {
                        Xi += 1.0;
                        S   = S*Xi/Lam + 1.0;
                    }
                    if (S < 0.0) X -= 1.0;
                }
            }
            return X;
        }
    }

    // bottom-up summation

    double Del;

    X   = 0.0;
    T   = exp(0.5*Lam);
    Del = 0.0;
    if (U > 0.5) Del = T*(1e-13*T);
    S   = 1.0 - T*(U*T) + Del;

    while (S < 0.0) {
        X  += 1.0;
        T   = X/Lam;
        Del = T*Del;
        S   = T*S + 1.0;
    }

    // top-down summation if needed

    if (S < 2.0*Del) {
        Del = 1e13*Del;
        T   = 1e17*Del;
        Del = (1.0-U)*Del;

        while (Del < T) {
            X   += 1.0;
            Del *= X/Lam;
        }

        S = Del;
        T = 1.0;
        while (S > 0.0) {
            T *= X/Lam;
            S -= T;
            X -= 1.0;
        }
    }
    return X;
}

// rlog1(x) = x - log(1+x) for -0.618<x<1.618 using rational minimax approx.
__device__ static inline double rlog1(double X) {

#define P0  0.2000000000000000
#define P1 -0.3636535967811319
#define P2  0.2015244511825799
#define P3 -0.03274937605228191
#define P4  0.00004542775258423288
#define Q1 -2.532553698191349
#define Q2  2.261033627633934
#define Q3 -0.8263420776370621
#define Q4  0.1008870710149048

    double R, T, W;

    // rlog1(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    R = X/(X + 2.0);
    T = R*R;
    W = (1.0/3.0) + T*((((P4*T + P3)*T + P2)*T + P1)*T + P0 ) /
                      ((((Q4*T + Q3)*T + Q2)*T + Q1)*T + 1.0);
    return T*((X + 2.0) - 2.0*R*W);
} /* rlog1 */
#undef P0
#undef P1
#undef P2
#undef P3
#undef P4
#undef Q1
#undef Q2
#undef Q3
#undef Q4
