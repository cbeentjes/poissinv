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

#ifndef POISSINV_H
#define POISSINV_H

// Standard Math Library or Intel Math Library
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

// Mathematical constants definition
#ifndef M_2PI
#define M_2PI 6.2831853071795865 /* 2*pi */
#endif

// declare prototype for inverse Normal CDF function
// defined at the bottom of this header file

static double normcdfinv_as241(double);
#define norminv(u) normcdfinv_as241(u);

// declare prototype for rlog1(x) = x - log(1+x) on -0.6<x<1.6
static double rlog1(double);

// As described in the TOMS paper, there are two versions;
// the first  (poissinv) is optimised for MIMD execution, whereas
// the second (poissinv_v) is designed for SIMD/vector execution

// poissinv(u, lam):
//
// This double precision function computes the inverse
// of the Poisson CDF
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~9e+15
// ave |error| < 1.6e-15*sqrt(max(1e-2,lam)) for 0    <= lam <= 13.5
//             < 3.0e-16*sqrt(lam)           for 13.5 <  lam <  ~9e+15
//
// For lam >= ~9e+15, the errors will be about 1 ulp.
//

static double poissinv_core(double, double, double);

double poissinv(double U, double Lam) {
    return poissinv_core(U, 1.0-U, Lam);
}

double poisscinv(double V, double Lam) {
    return poissinv_core(1.0-V, V, Lam);
}

// Forcing inline here improves performance typically by ~1% to 5%
inline __attribute__((always_inline)) static double poissinv_core(double U, double V, double Lam) {

    double X=0.0, Xi, W, T, Del, R, R2, S, S2, Eta, B0, B1, Lami=1.0/Lam;

    // handle exceptions -- constants defined in <math.h>

    // handles NAN inputs as well
    if (!(Lam >= 0.0) || isinf(Lam)) return FP_NAN;
    if (!(U >= 0.0 && V >= 0.0)) return FP_NAN;
    if (U == 0.0 || Lam == 0.0)  return 0.0;
    if (V == 0.0) return FP_INFINITE;

    if (Lam > 13.5) {

        double Lam_root = sqrt(Lam);

        W = norminv(fmin(U,V));
        if (U > V) W = -W;

// use polynomial approximations in central region

        if (fabs(W) < 3.0) {;

            // Q_N2
            S = Lam_root*W + (1.0/3.0 + (1.0/6.0)*W*W)*(1.0 - W/(12.0*Lam_root));

            Del = (1.0 /160.0);
            Del = (1.0 / 80.0) + Del*(W*W);
            Del = (1.0 / 40.0) + Del*(W*W);
            Del = Del * Lami;

            // Sum from smallest to largest to minimise rounding error;
            S += Del;
        }

// otherwise use Temme uniform approximation

        else {
            double  Rm=-1.0, C0=0.0;

            S = W / Lam_root;
            // use asymptotic expansions to solve for r when eta near 0
            if (fabs(S) < 1e-2) {
                // expansion of f^{-1}(s) - 1
                Rm =  (19.0/34020.0);
                Rm = -(23.0/17280.0) + Rm*S;
                Rm =  ( 1.0/  270.0) + Rm*S;
                Rm = -( 1.0/   72.0) + Rm*S;
                Rm =  ( 1.0/    6.0) + Rm*S;
                Rm =    1.0          + Rm*S;
                Rm *= S;
                R  = 1.0 + Rm;

                // O(1) correction
                C0 = -(137.0/38880.0);
                C0 =  (  7.0/  810.0) + C0*S;
                C0 = -(  1.0/   36.0) + C0*S;
                C0 =  (  1.0/    3.0) + C0*S;

                Del = 2.5e-2*Lami;
            }
            // use Newton iteration to solve for r when eta not near 0
            else if (S > -sqrt(2.0)) {
                R = 1.0 + S;
                if (R < 0.1) R = 0.1;

                do {
                    T  = log(R);
                    R2 = R;
                    S2 = sqrt(2.0*((1.0-R) + R*T));
                    if (R < 1.0) S2 = -S2;
                    R = R2 - (S2-S)*S2/T;
                    if (R < 0.1*R2) R = 0.1*R2;
                } while (fabs(R-R2) > 1e-6);

                T = log(R);
                Rm = R - 1.0;
                C0 = log(sqrt(2.0*R*((1.0-R)+R*T))/fabs(R-1.0)) / T;

                // O(1/lam) ad-hoc correction
                C0 += -0.0218*Lami/(0.065 + R);
                Del = 0.01*Lami/R;

                // If not using O(1/lam) correction delta is larger
                // Del = 0.03*Lami/R;
            }
            // Otherwise in the rare scenario where Temme approximation cannot
            // yield a solution. Make sure that direct summation is used
            // by setting S = -Lam so that X=trunc(S+Lam)=0.0
            else {
                Del= 0.0;
            }
            // Sum from smallest to largest to minimise rounding error;
            S = (Del + C0) + Lam*Rm;
        }

        // NOTE: use of round down mode in final sum is important for
        // large values of Lam to ensure fast execution
        // If Lam >> 1 then Del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        // Only apply the round down transformation if needed as this
        // operation is (relatively) costly on CPUs and not often needed.

        // Round down to nearest integer and check accuracy if X > 10

        if (Lam < 1e+15 || Lam > 1e+16) {
            X = trunc(S+Lam);
        } else {
            X = S;
            X += Lam;
            if ((X - Lam) >= S) {      // If S += Lam was rounded up
                X = nextafter(X,0.0);  // use round down mode instead
            }
            X = trunc(X);
        }

        if (X > 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0*Del >= X - Lam || X > 9.007199254740992e+15) return X;

            // Not accurate enough, need to check the CDF

            // correction procedure based on Temme approximation

            if (X > 0.5*Lam && X < 2.0*Lam) {

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

                S = S*exp(-0.5*X*Eta*Eta)/sqrt(M_2PI*X);
                if (X < Lam) {
                    S += 0.5*erfc(Eta*sqrt(0.5*X));
                    if (S > U) X -= 1.0;
                }
                else {
                    S -= 0.5*erfc(-Eta*sqrt(0.5*X));
                    if (S > -V) X -= 1.0;
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
                S = (X - Lam) - X*log(X*Lami) - S;

                // downwards

                if (X < Lam) {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*(U*T) * sqrt(M_2PI*Xi) * Lam;
                    T  = 1.0;
                    Xi = X;
                    for (int i=0; i < 50; i++) {
                        Xi -= 1.0;
                        T  *= Xi*Lami;
                        S  += T;
                    }
                    if (S > 0.0) X -= 1.0;
                }
                // upwards
                else {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*(V*T) * sqrt(M_2PI*X);
                    Xi = X;
                    for (int i=0; i<50; i++) {
                        Xi += 1.0;
                        S   = S*Xi*Lami + 1.0;
                    }
                    if (S < 0.0) X -= 1.0;
                }
            }
            return X;
        }
    }

    // bottom-up summation

    X   = 0.0;
    T   = exp(0.5*Lam);
    Del = 0.0;
    if (U > 0.5) Del = T*(1e-13*T);
    S   = 1.0 - T*(U*T) + Del;

    while (S < 0.0) {
        X  += 1.0;
        T   = X*Lami;
        Del = T*Del;
        S   = T*S + 1.0;
    }

    // top-down summation if needed

    if (S < 2.0*Del) {
        Del = 1e13*Del;
        T   = 1e17*Del;
        Del = V*Del;

        while (Del < T) {
            X   += 1.0;
            Del *= X*Lami;
        }

        S = Del;
        T = 1.0;
        while (S > 0.0) {
            T *= X*Lami;
            S -= T;
            X -= 1.0;
        }
    }
    return X;
}

// poissinv_v(u, lam):
//
// This double precision function computes the inverse
// of the Poisson CDF using an approach more suitable for SIMD execution
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~9e+15
// ave |error| < 4.5e-16*max(1,lam)   for 0    <= lam < 13.5
//             < 2.6e-16*sqrt(lam)    for 13.5 <  lam < ~9e+15
//
// For lam >= ~9e+15, the errors will be about 1 ulp.
//

double poissinv_v(double U, double Lam) {

    double X=0.0, Xi, T, Del, Rm, R, R2, S, S2, Eta, B0, B1;

    // handle exceptions -- constants defined in <math.h>

    if (!(Lam >= 0.0) || isinf(Lam)) return FP_NAN;
    if (!(U >= 0.0 && U <= 1.0)) return FP_NAN;
    if (U == 0.0 || Lam == 0.0)  return 0.0;
    if (U == 1.0) return FP_INFINITE;

    // large lam

    if (Lam > 4.0) {
        double W;

        W = norminv(U);

        S   = W/sqrt(Lam);
        Del = 3.3e-6;

        // use polynomial approximations in central region

        if ((S > -0.6833501) && (S < 1.777993)) {;

            //  polynomial approximation to f^{-1}(s) - 1

            Rm =    2.82298751e-07;
            Rm =   -2.58136133e-06 + Rm*S;
            Rm =    1.02118025e-05 + Rm*S;
            Rm =   -2.37996199e-05 + Rm*S;
            Rm =    4.05347462e-05 + Rm*S;
            Rm =   -6.63730967e-05 + Rm*S;
            Rm =    0.000124762566 + Rm*S;
            Rm =   -0.000256970731 + Rm*S;
            Rm =    0.000558953132 + Rm*S;
            Rm =    -0.00133129194 + Rm*S;
            Rm =     0.00370367937 + Rm*S;
            Rm =     -0.0138888706 + Rm*S;
            Rm = 0.166666666666666 + Rm*S;

            S +=                  S*(Rm*S);
            Rm = S;

            //  polynomial approximation to correction c0(r)

            T  =   1.86386867e-05;
            T  =  -0.000207319499 + T*Rm;
            T  =     0.0009689451 + T*Rm;
            T  =   -0.00247340054 + T*Rm;
            T  =    0.00379952985 + T*Rm;
            T  =   -0.00386717047 + T*Rm;
            T  =    0.00346960934 + T*Rm;
            T  =   -0.00414125511 + T*Rm;
            T  =    0.00586752093 + T*Rm;
            T  =   -0.00838583787 + T*Rm;
            T  =     0.0132793933 + T*Rm;
            T  =     -0.027775536 + T*Rm;
            T  =      0.333333333 + T*Rm;

            //  O(1/lam) correction

            X  =  -0.000134549989;
            X  =    0.00137126732 + X*Rm;
            X  =    -0.0058318548 + X*Rm;
            X  =     0.0134526835 + X*Rm;
            X  =    -0.0185770659 + X*Rm;
            X  =     0.0169808982 + X*Rm;
            X  =    -0.0135302102 + X*Rm;
            X  =     0.0135636388 + X*Rm;
            X  =      -0.01551508 + X*Rm;
            X  =     0.0174051003 + X*Rm;
            X  =    -0.0198016502 + X*Rm;
            X  =  X / Lam;

            //    sum from smallest to largest to minimise rounding error

            S = ((X + Del) + T) + Lam*Rm;
        }

        // otherwise use Newton iteration

        else if (S > -sqrt(2.0)) {

            R = 1.0 + S;
            if (R < 0.1) R = 0.1;

            do {
                T  = log(R);
                R2 = R;
                S2 = sqrt(2.0*(1.0 - R + R*T));
                if (R < 1.0) S2 = -S2;
                R = R2 - (S2-S)*S2/T;
                if (R < 0.1*R2) R = 0.1*R2;
            } while (fabs(R-R2) > 1e-5);

            T   = log(R);
            Rm = R - 1.0;
            T = log(sqrt(2.0*R*((1.0-R)+R*T))/fabs(R-1.0)) / T;

            // O(1/lam) ad-hoc correction
            T += -0.0218/(Lam*(0.065 + R));
            Del = 0.01/(Lam*R);

            // If not using O(1/lam) correction delta is larger
            // Del = 0.03/(Lam*R);

            S = (Del + T) + Lam*Rm;
        }

        // NOTE: use of round down mode in final sum is important for
        // large values of lam to ensure fast execution
        // If Lam >> 1 then del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        // On GPU/SIMD devices simply always use this instruction.

        // Round down to nearest integer and check accuracy if X > 10

        // This emulates the use of round-down mode, rather than relying on
        // the specific instruction
        X = S + Lam;
        if ( (X - Lam) > S ) {
            X = nextafter(X,0.0);
        }
        X = trunc(X);

        if (X > 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0*Del >= X - Lam || X > 9.007199254740992e+15) return X;

            // Not accurate enough, need to check the CDF

            // correction procedure based on Temme approximation

            if (X > 0.5*Lam && X < 2.0*Lam) {

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

                S = S*exp(-0.5*X*Eta*Eta)/sqrt(M_2PI*X);
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
                    S  = 1.0 - T*(U*T) * sqrt(M_2PI*Xi) * Lam;
                    T  = 1.0;
                    Xi = X;
                    for (int i=0; i<50; i++) {
                        Xi -= 1.0;
                        T  *= Xi/Lam;
                        S  += T;
                    }
                    if (S > 0.0) X -= 1.0;
                }

                else {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*((1.0-U)*T) * sqrt(M_2PI*X);
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
static double rlog1(double x) {

#define P0  0.2000000000000000
#define P1 -0.3636535967811319
#define P2  0.2015244511825799
#define P3 -0.03274937605228191
#define P4  0.00004542775258423288
#define Q1 -2.532553698191349
#define Q2  2.261033627633934
#define Q3 -0.8263420776370621
#define Q4  0.1008870710149048

    double r, t, w;

    // rlog1(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    r = x/(x + 2.0);
    t = r*r;
    w = (1.0/3.0) + t*((((P4*t + P3)*t + P2)*t + P1)*t + P0 ) /
                      ((((Q4*t + Q3)*t + Q2)*t + Q1)*t + 1.0);
    return t*((x + 2.0) - 2.0*r*w);
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

//////////////////////////////////////////////////////////////////////
//                                                                  //
// The routine below is a C version of the code in                  //
//                                                                  //
// ALGORITHM AS241: APPLIED STATS (1988) VOL. 37, NO. 3, 477-44.    //
// http://lib.stat.cmu.edu/apstat/241                               //
//                                                                  //
// The relative error is less than 1e-15, and the accuracy is       //
// verified in the accompanying MATLAB code as241_test.m            //
//                                                                  //
//////////////////////////////////////////////////////////////////////

static double normcdfinv_as241(double p) {

    double q, r, num, den, res;

    q = p - 0.5;
    if (fabs(q) <= 0.425) {
        r = 0.180625 - q*q;

        num =         2.5090809287301226727e+3;
        num = r*num + 3.3430575583588128105e+4;
        num = r*num + 6.7265770927008700853e+4;
        num = r*num + 4.5921953931549871457e+4;
        num = r*num + 1.3731693765509461125e+4;
        num = r*num + 1.9715909503065514427e+3;
        num = r*num + 1.3314166789178437745e+2;
        num = r*num + 3.3871328727963666080e0;

        den =         5.2264952788528545610e+3;
        den = r*den + 2.8729085735721942674e+4;
        den = r*den + 3.9307895800092710610e+4;
        den = r*den + 2.1213794301586595867e+4;
        den = r*den + 5.3941960214247511077e+3;
        den = r*den + 6.8718700749205790830e+2;
        den = r*den + 4.2313330701600911252e+1;
        den = r*den + 1.0000000000e+00;

        res = q * num / den;

        return res;
    }

    else {

        if (q < 0.0)
            r = p;
        else
            r = 1.0 - p;

        r = sqrt(-log(r));

        if (r <= 5.0) {
            r = r - 1.6;

            num =         7.74545014278341407640e-4;
            num = r*num + 2.27238449892691845833e-2;
            num = r*num + 2.41780725177450611770e-1;
            num = r*num + 1.27045825245236838258e0;
            num = r*num + 3.64784832476320460504e0;
            num = r*num + 5.76949722146069140550e0;
            num = r*num + 4.63033784615654529590e0;
            num = r*num + 1.42343711074968357734e0;

            den =         1.05075007164441684324e-9;
            den = r*den + 5.47593808499534494600e-4;
            den = r*den + 1.51986665636164571966e-2;
            den = r*den + 1.48103976427480074590e-1;
            den = r*den + 6.89767334985100004550e-1;
            den = r*den + 1.67638483018380384940e0;
            den = r*den + 2.05319162663775882187e0;
            den = r*den + 1.0000000000e+00;

            res = num / den;
        }

        else {
            r = r - 5.0;

            num =         2.01033439929228813265e-7;
            num = r*num + 2.71155556874348757815e-5;
            num = r*num + 1.24266094738807843860e-3;
            num = r*num + 2.65321895265761230930e-2;
            num = r*num + 2.96560571828504891230e-1;
            num = r*num + 1.78482653991729133580e0;
            num = r*num + 5.46378491116411436990e0;
            num = r*num + 6.65790464350110377720e0;

            den =         2.04426310338993978564e-15;
            den = r*den + 1.42151175831644588870e-7;
            den = r*den + 1.84631831751005468180e-5;
            den = r*den + 7.86869131145613259100e-4;
            den = r*den + 1.48753612908506148525e-2;
            den = r*den + 1.36929880922735805310e-1;
            den = r*den + 5.99832206555887937690e-1;
            den = r*den + 1.0000000000e+00;

            res = num / den;
        }

        if (q < 0.0)
            res = - res;

        return res;
    }
}

#endif /* POISSINV_H */
