//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__); fprintf(fp, __VA_ARGS__);}

//
// poissinv header file
//

#include "poissinv.h"
#include "poissinvf.h"

//
// linux timing routine
//

#include <sys/time.h>

double elapsed_time(double *et) {
    struct timeval t;
    double old_time = *et;

    gettimeofday( &t, (struct timezone *)0 );
    *et = t.tv_sec + t.tv_usec*1.0e-6;

    return *et - old_time;
}

//
// function prototype for quad precision evaluation of poisson inverse CDF
//

void poissinv_quad(int, float, float*, float*, double*, double*);
void poisscinv_quad(int, float, float*, float*, double*, double*);

//
// functions to calculate the floating point jump points of the approximate
// inverse CDF of the poisson distribution
//


// Double precision CPU/MIMD code version
void poissinv_bisection_scalar(int n, int N, double lam,
                         double *ulo_d, double *uhi_d ) {

    double x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0;
        u_lo  = 0.0;
        u_mid = 0.5*(u_hi + u_lo);
        xt    = (double) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
        x = poissinv(u_mid, lam);

        if (x>xt)
          u_hi = u_mid;
        else
          u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

// Double precision GPU/SIMD code version
void poissinv_bisection_vector(int n, int N, double lam,
                         double *ulo_d, double *uhi_d ) {

    double x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0;
        u_lo  = 0.0;
        u_mid = 0.5*(u_hi + u_lo);
        xt    = (double) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
        x = poissinv_v(u_mid, lam);

        if (x>xt)
          u_hi = u_mid;
        else
          u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

// Single precision CPU/MIMD code version
void poissinvf_bisection_scalar(int n, int N, float lam,
                         float *ulo_d, float *uhi_d ) {

    float x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0f;
        u_lo  = 0.0f;
        u_mid = 0.5f*(u_hi + u_lo);
        xt    = (float) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
        x = poissinvf(u_mid, lam);

        if (x>xt)
          u_hi = u_mid;
        else
          u_lo = u_mid;

        u_mid = 0.5f*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

// Single precision GPU/SIMD code version
void poissinvf_bisection_vector(int n, int N, float lam,
                         float *ulo_d, float *uhi_d ) {

    float x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0f;
        u_lo  = 0.0f;
        u_mid = 0.5f*(u_hi + u_lo);
        xt    = (float) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
        x = poissinvf_v(u_mid, lam);

        if (x>xt)
          u_hi = u_mid;
        else
          u_lo = u_mid;

        u_mid = 0.5f*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

//
// functions to calculate the floating point jump points of the approximate
// inverse complimentary CDF of the poisson distribution
//

// Double precision CPU/MIMD code version
void poisscinv_bisection_scalar(int n, int N, double lam,
                         double *ulo_d, double *uhi_d ) {

    double x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0;
        u_lo  = 0.0;
        u_mid = 0.5*(u_hi + u_lo);
        xt    = (double) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
        x = poisscinv(u_mid, lam);

        if (x<xt)
          u_hi = u_mid;
        else
          u_lo = u_mid;

        u_mid = 0.5*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

// Single precision CPU/MIMD code version
void poisscinvf_bisection_scalar(int n, int N, float lam,
                         float *ulo_d, float *uhi_d ) {

    float x, xt, u_lo, u_hi, u_mid;

    if (n <= N) {
        u_hi  = 1.0f;
        u_lo  = 0.0f;
        u_mid = 0.5f*(u_hi + u_lo);
        xt    = (float) n;

        while ((u_mid>u_lo) & (u_mid<u_hi)) {
            x = poisscinvf(u_mid, lam);

            if (x<xt)
                u_hi = u_mid;
            else
                u_lo = u_mid;

            u_mid = 0.5f*(u_hi + u_lo);
        }
        ulo_d[n] = u_lo;
        uhi_d[n] = u_hi;
    }
}

//
// function to compare the poissinv floating point jump points to exact jump
// points calculated using quad precision
//

void poissinv_err_compare(FILE *fp, int N, float lam,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2;

    float lam_root = sqrtf(lam);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision
    //////////////////////////////////////////////////

    elapsed_time(&timer);
    poissinv_quad(N, lam, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1  = 0.0;
    err1f = 0.0;

    elapsed_time(&timer);
    for (int n=0; n<=N; n++) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by poissinv(f)
        if (n < lam - 200.0 - 50.0*lam_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        // Similarly, no need to continue if we are already beyond
        // double precision limit.
        else if (n > lam + 20.0 + 10.0*lam_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        }
        else {
            poissinv_bisection_scalar(n, N, (double) lam, Ulo_h, Uhi_h);
            poissinvf_bisection_scalar(n, N, lam, ulo_h, uhi_h);
        }
#else
        poissinv_bisection_scalar(n, N, (double) lam, Ulo_h, Uhi_h);
        poissinvf_bisection_scalar(n,(float) N, lam, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 16777216) {
            if ( uhi_h[n] <= ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, uhi_h[n] = %20.16g, ulo_ex[n-1] = %20.16g \n",
                        lam,      n,      uhi_h[n],           ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
            if ( uhi_ex[n] <= ulo_h[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                        lam,      n,      uhi_ex[n],           ulo_h[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                        lam,      n,      Uhi_h[n],           Ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
            if ( Uhi_ex[n] <= Ulo_h[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                        lam,      n,      Uhi_ex[n],           Ulo_h[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }

#ifdef FAST_CHECK
        // Check if we can break early
        if (n > lam + 20.0 + 20.0*lam_root) {
            break;
        }
#endif
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 error
    PRINTF2(fp,"%11.5f    %#9.2e     %#9.2e       %.2e    %.2e\n",lam,err1f,err1,elapsed1,elapsed2);
}

//
// function to compare the poissinv_v floating point jump points to exact jump
// points calculated using quad precision
//

void poissinv_v_err_compare(FILE *fp, int N, float lam,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2;

    float lam_root = sqrtf(lam);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision
    //////////////////////////////////////////////////

    elapsed_time(&timer);
    poissinv_quad(N, lam, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1  = 0.0;
    err1f = 0.0;

    elapsed_time(&timer);
    for (int n=0; n<=N; n++) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by poissinv(f)_v
        if (n < lam - 200.0 - 50.0*lam_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        else if (n > lam + 20.0 + 10.0*lam_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        }
        // Similarly, no need to continue if we are already beyond
        // double precision limit.
        else {
            poissinv_bisection_vector(n, N, (double) lam, Ulo_h, Uhi_h);
            poissinvf_bisection_vector(n, N, lam, ulo_h, uhi_h);
        }
#else
        poissinv_bisection_vector(n, N, (double) lam, Ulo_h, Uhi_h);
        poissinvf_bisection_vector(n,(float) N, lam, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 16777216) {
            if ( uhi_h[n] <= ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, uhi_h[n] = %20.16g, ulo_ex[n-1] = %20.16g \n",
                        lam,      n,      uhi_h[n],           ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
            if ( uhi_ex[n] <= ulo_h[n-1] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                        lam,      n,      uhi_ex[n],           ulo_h[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                        lam,      n,      Uhi_h[n],           Ulo_ex[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
            if ( Uhi_ex[n] <= Ulo_h[n-1] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                        lam,      n,      Uhi_ex[n],           Ulo_h[n-1]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }

#ifdef FAST_CHECK
        // Check if we can break early
        if (n > lam + 20.0 + 20.0*lam_root) {
            break;
        }
#endif
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 error
    PRINTF2(fp,"%11.5f    %#9.2e     %#9.2e       %.2e    %.2e\n",lam,err1f,err1,elapsed1,elapsed2);
}

//
// function to compare the poisscinv floating point jump points to exact jump
// points calculated using quad precision
//

void poisscinv_err_compare(FILE *fp, int N, float lam,
        double *Ulo_ex, double *Uhi_ex, double *Ulo_h, double *Uhi_h,
        float  *ulo_ex, float  *uhi_ex, float  *ulo_h, float  *uhi_h) {

    double err1, err1f;
    double timer, elapsed1, elapsed2;

    float lam_root = sqrtf(lam);

    //////////////////////////////////////////////////
    // compute reference solution in quad precision
    //////////////////////////////////////////////////

    elapsed_time(&timer);
    poisscinv_quad(N, lam, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);
    elapsed1 = elapsed_time(&timer);

    //////////////////////////////////////////////////
    //  check double and single precision versions
    //////////////////////////////////////////////////

    err1  = 0.0;
    err1f = 0.0;

    elapsed_time(&timer);
    for (int n=0; n<=N; n++) {
#ifdef FAST_CHECK // Speed up check by using knowledge of distribution shape
        // No need to continue if we are already beyond
        // double precision limit.
        if (n > lam + 200.0 + 50.0*lam_root) {
            Ulo_h[n] = 0.0;
            Uhi_h[n] = DBL_TRUE_MIN;
            ulo_h[n] = 0.0f;
            uhi_h[n] = FLT_TRUE_MIN;
        }
        // Determine good starting point of the summation, as we do not need
        // to evaluate for n so small that it cannot be generated by poisscinv(f)
        else if (n < lam - 10.0 - 10.0*lam_root){
            Ulo_h[n] = 1.0-DBL_EPSILON/2.0;
            Uhi_h[n] = 1.0;
            ulo_h[n] = 1.0f-FLT_EPSILON/2.0f;
            uhi_h[n] = 1.0f;
        }
        else {
            poisscinv_bisection_scalar(n, N, (double) lam, Ulo_h, Uhi_h);
            poisscinvf_bisection_scalar(n, N, lam, ulo_h, uhi_h);
        }
#else
        poisscinv_bisection_scalar(n, N, (double) lam, Ulo_h, Uhi_h);
        poisscinvf_bisection_scalar(n, N, lam, ulo_h, uhi_h);
#endif /* FAST_CHECK */

        // Single precision
        err1f += 0.5f*fabsf( (ulo_ex[n]-ulo_h[n]) + (uhi_ex[n]-uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 16777216) {
            if ( uhi_h[n-1] <= ulo_ex[n] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, uhi_h[n-1] = %20.16g, ulo_ex[n] = %20.16g \n",
                        lam,      n,      uhi_h[n-1],           ulo_ex[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
            if ( uhi_ex[n-1] <= ulo_h[n] ) {
                PRINTF2(fp,"\n Single precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, uhi_ex[n-1] = %20.16g, ulo_h[n] = %20.16g \n",
                        lam,      n,      uhi_ex[n-1],           ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",ulo_h[n], uhi_h[n], ulo_ex[n-1], uhi_ex[n-1]);
                exit(1);
            }
        }

        // Double precision
        err1  += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

        // check if Linfinity error is <= 1. Note that this is only relevant
        // if ULP < 1
        if (n > 0 && n < 9007199254740992) {
            if ( Uhi_h[n-1] <= Ulo_ex[n] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_h[n-1] = %20.16g, Ulo_ex[n] = %20.16g \n",
                        lam,      n,      Uhi_h[n-1],           Ulo_ex[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
            if ( Uhi_ex[n-1] <= Ulo_h[n] ) {
                PRINTF2(fp,"\n Double precision error \n");
                PRINTF2(fp,"\n error: lam = %f, n = %d, Uhi_ex[n-1] = %20.16g, Ulo_h[n] = %20.16g \n",
                        lam,      n,      Uhi_ex[n-1],           Ulo_h[n]);
                PRINTF2(fp,"%20.16g, %20.16g, %20.16g, %20.16g\n",Ulo_h[n], Uhi_h[n], Ulo_ex[n-1], Uhi_ex[n-1]);
                exit(1);
            }
        }

#ifdef FAST_CHECK
        // Check if we can break early
        if (n > lam + 200.0 + 150.0*lam_root) {
            break;
        }
#endif
    }
    elapsed2 = elapsed_time(&timer);
    // print parameters and L1 error
    PRINTF2(fp,"%11.5f    %#9.2e     %#9.2e       %.2e    %.2e\n",lam,err1f,err1,elapsed1,elapsed2);
}


//////////////////////////////////////////////////
// main code
//////////////////////////////////////////////////

int main(int argc, char **argv) {

    // Destination file
    char filename[128];
    FILE *fp;
    int N;

    if (argc > 2) {
        sprintf(filename, "%s", argv[2]);
    } else {
        sprintf(filename, "poissinv_check.txt");
    }
    fp = fopen(filename, "w");

    // Single precision variables
    float  *ulo_h, *uhi_h, *ulo_ex, *uhi_ex;
    // Double precision variables
    double *Ulo_h, *Uhi_h, *Ulo_ex, *Uhi_ex;


    // Poisson distribution parameters
    float lam;

    float dl;
    float lam_max;
#ifdef PLUS
    dl = ldexpf(1, -13);
    lam_max = 1.1e+6f;
#else
    dl = 2.00f;
    lam_max = 2.0e+7f;
#endif

    // allocate memory
    int Nmax = (1<<30);

    ulo_ex = (float  *)malloc(Nmax*sizeof(float));
    uhi_ex = (float  *)malloc(Nmax*sizeof(float));
    ulo_h  = (float  *)malloc(Nmax*sizeof(float));
    uhi_h  = (float  *)malloc(Nmax*sizeof(float));
    Ulo_ex = (double *)malloc(Nmax*sizeof(double));
    Uhi_ex = (double *)malloc(Nmax*sizeof(double));
    Ulo_h  = (double *)malloc(Nmax*sizeof(double));
    Uhi_h  = (double *)malloc(Nmax*sizeof(double));

    // number of loop iterations in accuracy tests
#ifndef COUNT_LAMBDA
#define COUNT_LAMBDA 16
#endif

#ifndef SKIP_SCALAR

    // set values to test
    lam = 0.5f;
    // change lambda based on input if desired
    if (argc > 1) {
        lam = strtod(argv[1], NULL);
    }
    N   = 50 + (int) (2.0f*lam);
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests \n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"     lambda    error(FP32)   error(FP64)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_LAMBDA; count++) {
        if (N >= Nmax) exit(1);
        poissinv_err_compare(fp, N, lam, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        lam += dl;
#else
        lam *= dl;
#endif
        if (lam > lam_max) break;

        N   = 50 + (int) (2*lam);
    }
#endif /* SKIP_SCALAR */

    //////////////////////////////////////////////////
    // re-do for vector version
    //////////////////////////////////////////////////

#ifdef VECTOR

    // set values to test
    lam = 0.5f;
    // change lambda based on input if desired
    if (argc > 1) {
        lam = strtod(argv[1], NULL);
    }
    N   = 50 + (int) (2.0f*lam);
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nvector CPU algorithm accuracy tests \n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"     lambda    error(FP32)   error(FP64)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_LAMBDA; count++) {
        if (N >= Nmax) exit(1);
        poissinv_v_err_compare(fp, N, lam, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        lam += dl;
#else
        lam *= dl;
#endif
        if (lam > lam_max) break;

        N   = 50 + (int) (2*lam);
    }

#endif /* VECTOR */

    //////////////////////////////////////////////////
    // re-do for complementary version
    //////////////////////////////////////////////////

#ifdef COMPLEMENTARY

    PRINTF2(fp,"\n\n******************************\n");
    PRINTF2(fp,"Complementary version poisscinv\n");
    PRINTF2(fp,"******************************\n");

    // set values to test
    lam = 0.5f;
    // change lambda based on input if desired
    if (argc > 1) {
        lam = strtod(argv[1], NULL);
    }
    N   = 50 + (int) (2*lam);
    PRINTF2(fp,"\n");
    PRINTF2(fp,"\nscalar CPU algorithm accuracy tests \n");
    PRINTF2(fp,"----------------------------------------------------------------\n");
    PRINTF2(fp,"     lambda    error(FP32)   error(FP64)    T_quad(s)   T_num(s) \n");
    for (int count=0; count<COUNT_LAMBDA; count++) {
        poisscinv_err_compare(fp, N, lam, Ulo_ex, Uhi_ex, Ulo_h, Uhi_h, ulo_ex, uhi_ex, ulo_h, uhi_h);
#ifdef PLUS
        lam += dl;
#else
        lam *= dl;
#endif
        if (lam > lam_max) break;

        N   = 1000 + (int) (2*lam);
        if (N >= Nmax) exit(1);
    }

#endif /* COMPLEMENTARY */

    // free memory

    free(Ulo_h);
    free(Uhi_h);
    free(ulo_h);
    free(uhi_h);

    free(ulo_ex);
    free(uhi_ex);
    free(Ulo_ex);
    free(Uhi_ex);

    fclose(fp);
}
