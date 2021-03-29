#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <quadmath.h>

#include <float.h>

//////////////////////////////////////////////////
// compute Poisson p.m.f in quad precision
//////////////////////////////////////////////////

__float128 poisson_pmf_quad(int n, float lam) {
    __float128 q_lam, q_n;

    q_lam = (__float128) lam;
    q_n   = (__float128) n;

    if (lam==0.0f) return( (n==0) ? 1.0q : 0.0q);
    if (n==0)      return expq( -q_lam);

    return expq( -q_lam + q_n*logq(q_lam) - lgammaq(q_n+1.0q) );
}

//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

void poissinv_quad(int N, float lam, float  *u_lo, float  *u_hi,
        double *U_lo, double *U_hi) {
    __float128 q_s;

    q_s   = 0.0q;

    double lam_root = sqrt(lam);

    for (int n=0; n<=N; n++) {

        // For large lambda use normal approximation to determine good
        // starting point of the summation, as we do not need to evaluate
        // the PMF when it is smaller than quad precision.
        if (lam > 25000. && n < lam - 152.0*lam_root) {
            u_lo[n] = 0.0f;
            u_hi[n] = FLT_TRUE_MIN;
            U_lo[n] = 0.0;
            U_hi[n] = DBL_TRUE_MIN;
            continue;
        }

        // Similarly, no need to continue summing if we are already beyond
        // double precision.
        if (q_s > 1.0 - DBL_EPSILON/2.0) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
            continue;
        }

        q_s += poisson_pmf_quad(n, lam);

        //
        // single precision interval bisection
        //
        if (q_s < FLT_TRUE_MIN ) {
            u_hi[n] = FLT_TRUE_MIN;
            u_lo[n] = 0.0f;
        } else if (q_s > 1.0f - FLT_EPSILON/2.0f) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
        } else {
            u_hi[n] = 1.0f;
            u_lo[n] = 0.0f;
            float u_mid = 0.5f*(u_hi[n] + u_lo[n]);

            while ( (u_mid>u_lo[n]) & (u_mid<u_hi[n]) ) {
                if ((__float128) u_mid > q_s)
                    u_hi[n] = u_mid;
                else
                    u_lo[n] = u_mid;

                u_mid = 0.5f*(u_hi[n] + u_lo[n]);
            }
        }

        //
        // double precision interval bisection
        //
        if (q_s < DBL_TRUE_MIN ) {
            U_hi[n] = DBL_TRUE_MIN;
            U_lo[n] = 0.0;
        } else if (q_s > 1.0 - DBL_EPSILON/2.0) {
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
        } else {
            U_hi[n] = 1.0;
            U_lo[n] = 0.0;
            double U_mid = 0.5*(U_hi[n] + U_lo[n]);

            while ( (U_mid>U_lo[n]) & (U_mid<U_hi[n]) ) {
                if ((__float128) U_mid > q_s)
                    U_hi[n] = U_mid;
                else U_lo[n] = U_mid;
                U_mid = 0.5*(U_hi[n] + U_lo[n]);
            }
        }
    }
}

void poisscinv_quad(int N, float lam, float  *u_lo, float  *u_hi,
        double *U_lo, double *U_hi) {
    __float128 q_s;

    q_s   = 0.0q;

    double lam_root = sqrt(lam);

    for (int n=N; n>=0; n--) {

        // For large lambda use normal approximation to determine good
        // starting point of the summation, as we do not need to evaluate
        // the PMF when it is smaller than quad precision.
        if (lam > 25000. && n > lam + 152.0*lam_root) {
            u_lo[n] = 0.0f;
            u_hi[n] = FLT_TRUE_MIN;
            U_lo[n] = 0.0;
            U_hi[n] = DBL_TRUE_MIN;
            continue;
        }

        // Similarly, no need to continue summing if we are already beyond
        // double precision.
        if (q_s > 1.0 - DBL_EPSILON/2.0) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
            continue;
        }

        q_s += poisson_pmf_quad(n, lam);

        //
        // single precision interval bisection
        //
        if (q_s < FLT_TRUE_MIN ) {
            u_hi[n] = FLT_TRUE_MIN;
            u_lo[n] = 0.0f;
        } else if (q_s > 1.0f - FLT_EPSILON/2.0f) {
            u_hi[n] = 1.0f;
            u_lo[n] = 1.0f-FLT_EPSILON/2.0f;
        } else {
            u_hi[n] = 1.0f;
            u_lo[n] = 0.0f;
            float u_mid = 0.5f*(u_hi[n] + u_lo[n]);

            while ( (u_mid>u_lo[n]) & (u_mid<u_hi[n]) ) {
                if ((__float128) u_mid > q_s)
                    u_hi[n] = u_mid;
                else
                    u_lo[n] = u_mid;

                u_mid = 0.5f*(u_hi[n] + u_lo[n]);
            }
        }

        //
        // double precision interval bisection
        //
        if (q_s < DBL_TRUE_MIN ) {
            U_hi[n] = DBL_TRUE_MIN;
            U_lo[n] = 0.0;
        } else if (q_s > 1.0 - DBL_EPSILON/2.0) {
            U_hi[n] = 1.0;
            U_lo[n] = 1.0-DBL_EPSILON/2.0;
        } else {
            U_hi[n] = 1.0;
            U_lo[n] = 0.0;
            double U_mid = 0.5*(U_hi[n] + U_lo[n]);

            while ( (U_mid>U_lo[n]) & (U_mid<U_hi[n]) ) {
                if ((__float128) U_mid > q_s)
                    U_hi[n] = U_mid;
                else
                    U_lo[n] = U_mid;

                U_mid = 0.5*(U_hi[n] + U_lo[n]);
            }
        }
    }
}
