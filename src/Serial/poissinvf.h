/////////////////////////////////////////////////////////////////
//                                                             //
// This software was written by Mike Giles, 2014               //
//                                                             //
// Minor tweaks by Casper Beentjes, 2021                       //
//                                                             //
// It is copyright University of Oxford, and provided under    //
// the terms of the GNU GPLv3 license:                         //
// http://www.gnu.org/licenses/gpl.html                        //
//                                                             //
// Commercial users wanting to use the software under a more   //
// permissive license, such as BSD, should contact the author: //
// mike.giles@maths.ox.ac.uk                                   //
//                                                             //
/////////////////////////////////////////////////////////////////

#ifndef POISSINVF_H
#define POISSINVF_H

// Standard Math Library or Intel Math library
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// Mathematical constants definition
#ifndef M_2PIF
#define M_2PIF 6.2831855f /* 2*pi */
#endif

// declare prototype for inverse Normal CDF function
// defined at the bottom of this header file

static float normcdfinvf_as241(float);
#define norminvf(u) normcdfinvf_as241(u);

// declare prototype for rlog1(x) = x - log(1+x) on -0.6<x<1.6
static float rlog1f(float);

// As described in the TOMS paper, there are two versions;
// the first  (poissinvf) is optimised for MIMD execution, whereas
// the second (poissinvf_v) is designed for SIMD/vector execution

// poissinvf(u, lam):
//
// This single precision function computes the inverse
// of the Poisson CDF
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~1.6e+07
// ave |error| < 6.8e-07*sqrt(max(1e-02,lam)   for 0    <= lam <= 11.5
//             < 1.8e-07*sqrt(lam)             for 11.5 <  lam <  ~1.6e+07
//
// For lam > ~1.6e+07, the errors will be about 1 ulp.
//

static float poissinvf_core(float, float, float);

float poissinvf(float u, float lam) {
    return poissinvf_core(u, 1.0f-u, lam);
}

float poisscinvf(float v, float lam) {
    return poissinvf_core(1.0f-v, v, lam);
}

// Forcing inline here improves performance typically by ~1% to 5%
inline __attribute__((always_inline)) static float poissinvf_core(float u, float v, float lam) {

    float x=0.0f, xi, w, t, del, r, r2, s, s2, eta, b0, b1, lami=1.0f/lam;

    // handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(lam >= 0.0f) || isinf(lam)) return FP_NAN;
    if (!(u >= 0.0f && v >= 0.0f)) return FP_NAN;
    if (u == 0.0f || lam == 0.0f)  return 0.0f;
    if (v == 0.0f) return FP_INFINITE;

    if (lam > 11.5f) {

        float lam_root = sqrtf(lam);

        w = norminvf(fminf(u,v));
        if (u > v) w = -w;

// use polynomial approximation in central region

        if (fabsf(w) < 3.0f) {;

            s = lam_root*w + (1.0f/3.0f + (1.0f/6.0f)*w*w)*(1.0f - w/(12.0f*lam_root));

            del = (1.0f /160.0f);
            del = (1.0f / 80.0f) + del*(w*w);
            del = (1.0f / 40.0f) + del*(w*w);
            del = del * lami;

            s += del;
        }

// otherwise use Temme approximation

        else {
            float rm, c0;

            s = w / lam_root;

            // use asymptotic expansions to solve for r when eta near 0
            if (fabsf(s) < 1e-2f) {
                // expansion of f^{-1}(s) - 1
                rm =  (1.0f/6.0f);
                rm =   1.0f       + rm*s;
                rm *= s;
                r  = 1.0f + rm;

                // O(1) correction
                c0 =  1.0f/3.0f;

                // expansion is accurate enough to have relative error
                // within machine epsilon.
                del = 3.0e-08f*lam*r;
            }
            // use Newton iteration to solve for r when eta not near 0
            else {
                r = 1.0f + s;
                if (r < 0.1f) r = 0.1f;

                do {
                    t  = logf(r);
                    r2 = r;
                    s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                    if (r < 1.0f) s2 = -s2;
                    r = r2 - (s2-s)*s2/t;
                    if (r < 0.1f*r2) r = 0.1f*r2;
                } while (fabsf(r-r2) > 1e-5f);

                t = logf(r);
                rm = r - 1.0f;
                c0 = logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;

                // O(1/lam) ad-hoc correction
                c0 += -0.0218f*lami/(0.065f + r);
                del = 0.002f*lami/r;

                // If not using O(1/lam) correction delta is larger
                // del = 0.03f*lami/r;
            }
            // Sum from smallest to largest to minimise rounding error;
            s = (del + c0) + lam*rm;
        }

        // NOTE: use of round down mode in final sum is important for
        // large values of lam to ensure fast execution
        // If lam >> 1 then del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        // Only apply the round down transformation if needed as this
        // operation is (relatively) costly on CPUs and not often needed.

        // Round down to nearest integer and check accuracy if x > 10

        if (lam < 1.5e+6f || lam > 2e+7f) {
            x = truncf(s+lam);
        } else {
            x = s;
            x += lam;
            if ((x - lam) >= s) {        // If x += lam was rounded up
                x = nextafterf(x,0.0f);  // use round down mode instead
            }
            x = truncf(x);
        }

        if (x > 10.0f) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return x;
#endif

            // Accurate enough, no need to check the CDF
            if (s - 2.0f*del >= x - lam || x > 1.6777216e7f) return x;

            // Not accurate enough, need to check the CDF

            // correction procedure based on Temme approximation

            if (x > 0.5f*lam && x < 2.0f*lam) {

                xi = 1.0f / x;
                // Use rlog1(z) = z - log(1+z) function, a fast rational
                // minimax approximation for -0.6 < z < 1.6. This is a lot
                // more accurate when x is near 0 (here x near lambda).
                eta = sqrtf(2.0f*rlog1f((lam-x)*xi));
                if (x > lam) eta = -eta;

                b1 = -4.3820360184533521e-09f;               s = b1;
                b0 =  1.0261809784240299e-08f;               s = b0 + s*eta;
                b1 =  6.7078535434015332e-09f + 13.0f*b1*xi; s = b1 + s*eta;
                b0 = -1.7665952736826086e-07f + 12.0f*b0*xi; s = b0 + s*eta;
                b1 =  8.2967113409530833e-07f + 11.0f*b1*xi; s = b1 + s*eta;
                b0 = -1.8540622107151585e-06f + 10.0f*b0*xi; s = b0 + s*eta;
                b1 = -2.1854485106799979e-06f +  9.0f*b1*xi; s = b1 + s*eta;
                b0 =  3.9192631785224383e-05f +  8.0f*b0*xi; s = b0 + s*eta;
                b1 = -0.00017875514403292177f +  7.0f*b1*xi; s = b1 + s*eta;
                b0 =  0.00035273368606701921f +  6.0f*b0*xi; s = b0 + s*eta;
                b1 =   0.0011574074074074078f +  5.0f*b1*xi; s = b1 + s*eta;
                b0 =   -0.014814814814814815f +  4.0f*b0*xi; s = b0 + s*eta;
                b1 =    0.083333333333333329f +  3.0f*b1*xi; s = b1 + s*eta;
                b0 =    -0.33333333333333331f +  2.0f*b0*xi; s = b0 + s*eta;
                s  = s / (1.0f + b1*xi);

                s = s*expf(-0.5f*x*eta*eta)/sqrtf(M_2PIF*x);
                if (x < lam) {
                    s += 0.5f*erfcf(eta*sqrtf(0.5f*x));
                    if (s > u) x -= 1.0f;
                }
                else {
                    s -= 0.5f*erfcf(-eta*sqrtf(0.5f*x));
                    if (s > -v) x -= 1.0f;
                }
            }

            // sum downwards or upwards

            else {
                xi = 1.0f / x;
                s = - (1.0f/360.0f);
                s =   (1.0f/12.0f)   + s*xi*xi;
                s =                    s*xi;
                s = (x - lam) - x*logf(x*lami) - s;

                // downwards
                if (x < lam) {
                    t  = expf(-0.5f*s);
                    s  = 1.0f - t*(u*t) * sqrtf(M_2PIF*xi) * lam;
                    t  = 1.0f;
                    xi = x;
                    for (int i=0; i<20; i++) {
                        xi -= 1.0f;
                        t  *= xi*lami;
                        s  += t;
                    }
                    if (s > 0.0f) x -= 1.0f;
                }
                // upwards
                else {
                    t  = expf(-0.5f*s);
                    s  = 1.0f - t*(v*t) * sqrtf(M_2PIF*x);
                    xi = x;
                    for (int i=0; i<20; i++) {
                        xi += 1.0f;
                        s   = s*xi*lami + 1.0f;
                    }
                    if (s < 0.0f) x -= 1.0f;
                }
            }
            return x;
        }
    }


    // bottom-up summation

    x   = 0.0f;
    t   = expf(0.5f*lam);
    del = 0.0f;
    if (u > 0.5f) del = t*(1e-6f*t);
    s   = 1.0f - t*(u*t) + del;

    while (s < 0.0f) {
        x  += 1.0f;
        t   = x*lami;
        del = t*del;
        s   = t*s + 1.0f;
    }

    // top-down summation if needed

    if (s < 2.0f*del) {
        del = 1e6f*del;
        t   = 1e7f*del;
        del = v*del;

        while (del < t) {
            x   += 1.0f;
            del *= x*lami;
        }

        s = del;
        t = 1.0f;
        while (s > 0.0f) {
            t *= x*lami;
            s -= t;
            x -= 1.0f;
        }
    }

    return x;
}

// poissinvf_v(u, lam):
//
// This single precision function computes the inverse
// of the Poisson CDF using an approach more suitable for SIMD
// execution
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~1.6e+07
// ave |error| < 6.2e-07*sqrt(max(1e-02,lam)   for 0  <= lam <= 6
//             < 1.8e-07*sqrt(max(69   ,lam)   for 6  <  lam <  ~1.6e+07
//
// For lam > ~1.6e+07, the errors will be about 1 ulp.
//

float poissinvf_v(float u, float lam) {

    float s, t, x = 0.0f;

    // handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(lam >= 0.0f) || isinf(lam)) return FP_NAN;
    if (!(u >= 0.0f && u <= 1.0f)) return FP_NAN;
    if (u == 0.0f || lam == 0.0f)  return 0.0f;
    if (u == 1.0f) return FP_INFINITE;

    if (lam > 4.0f) {

        float w;
        w = norminvf(u);
        s = w/sqrtf(lam);

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
            x  /= lam;

            // NOTE: use of round down mode in final sum is important for
            // large values of lambda to ensure correct rounding
            // On GPU/SIMD devices simply always use this instruction.

            s = (x+t) + lam*rm;
            x = s + lam;
            // Make sure to round down
            if ( (x - lam) > s) {
                x = nextafterf(x,0.0f);
            }
            x = truncf(x);
        }

        // otherwise use Newton iteration

        else if (s > -sqrtf(2.0f)) {
            float r, r2, s2;

            r = 1.0f + s;
            if (r < 0.1f) r = 0.1f;

            do {
                t  = logf(r);
                r2 = r;
                s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                if (r < 1.0f) s2 = -s2;
                r = r2 - (s2-s)*s2/t;
                if (r < 0.1f*r2) r = 0.1f*r2;
            } while (fabsf(r-r2) > 1e-5f);

            t = logf(r);
            x = lam*r + logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;
            x = truncf( x - 0.0218f/(x+0.065f*lam) );
        }
    }

    // bottom-up summation

    if (x < 10.0f) {
        float del;
        x   = 0.0f;
        t   = expf(0.5f*lam);
        del = 0.0f;
        if (u > 0.5f) del = t*(1e-6f*t);
        s   = 1.0f - t*(u*t) + del;

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

// rlog1f(x) = x - log(1+x) for -0.618<x<1.618 using rational minimax approx.
inline static float rlog1f(float x) {

#define P0  0.19999999f
#define P1 -0.11665086f
#define Q1 -1.2975472f
#define Q2  0.37142160f

    float r, t, w;

    // rlog1f(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    r = x/(x + 2.0f);
    t = r*r;
    w = (1.0f/3.0f) + t*(P1*t + P0) /
        ((Q2*t + Q1)*t + 1.0f);
    return t*((x + 2.0f) - 2.0f*r*w);
} /* rlog1f */
#undef P0
#undef P1
#undef Q1
#undef Q2

//////////////////////////////////////////////////////////////////////
//                                                                  //
// The routine below is a C version of the code in                  //
//                                                                  //
// ALGORITHM AS241: APPLIED STATS (1988) VOL. 37, NO. 3, 477-44.    //
// http://lib.stat.cmu.edu/apstat/241                               //
//                                                                  //
// The relative error is less than 5 ULP in single precision, and   //
// the accuracy is verified in the accompanying MATLAB code         //
// as241_test.m                                                     //
//                                                                  //
//////////////////////////////////////////////////////////////////////

static float normcdfinvf_as241(float p) {

    float q, r, num, den, res;

    q = p - 0.5f;
    if (fabsf(q) <= 0.425f) {
        r = 0.180625f - q*q;

        num =         5.9109374720e+1f;
        num = r*num + 1.5929113202e+2f;
        num = r*num + 5.0434271938e+1f;
        num = r*num + 3.3871327179e+0f;

        den =         6.7187563600e+1f;
        den = r*den + 7.8757757664e+1f;
        den = r*den + 1.7895169469e+1f;
        den = r*den + 1.0000000000e+0f;

        res = q * num / den;

        return res;
    }

    else {

        if (q < 0.0f)
            r = p;
        else
            r = 1.0f - p;

        r = sqrtf(-logf(r));

        if (r <= 5.0f) {
            r = r - 1.6f;

            num =         1.7023821103e-1f;
            num = r*num + 1.3067284816e+0f;
            num = r*num + 2.7568153900e+0f;
            num = r*num + 1.4234372777e+0f;

            den =         1.2021132975e-1f;
            den = r*den + 7.3700164250e-1f;
            den = r*den + 1.0000000000e+0f;

            res = num / den;
        }

        else {
            r = r - 5.0f;

            num =         1.7337203997e-2f;
            num = r*num + 4.2868294337e-1f;
            num = r*num + 3.0812263860e+0f;
            num = r*num + 6.6579051150e+0f;

            den =         1.2258502635e-2f;
            den = r*den + 2.4197894225e-1f;
            den = r*den + 1.0000000000e+0f;

            res = num / den;
        }

        if (q < 0.0f)
            res = - res;

        return res;
    }
}

#endif /* POISSINVF_H */
