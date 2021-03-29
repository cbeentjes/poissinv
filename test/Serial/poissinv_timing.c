//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__); fprintf(fp,  __VA_ARGS__);}

//
// poissinv header files
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
// CDF inverse timing test functions
//

void normcdfinv_test_double_scalar(int N) {
    double x, u;
    int n;

    for (n=0; n<N; n++) {
        // Generate "random" variates
        u = (n+0.5) / N;

        x = norminv(u);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Double precision CPU/MIMD-style version
void poissinv_test_double_scalar(int N, double lam) {
    double x, u;
    int    n;

    for (n=0; n<N; n++) {
        // Generate "random" variates
        u = (n+0.5) / N;

        // extra n*1e-100 is to prevent compiler
        // optimising for fixed lam
        x = poissinv(u, lam+n*1e-100);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Double precision GPU/SIMD-style version
void poissinv_test_double_vector(int N, double lam) {
    double x, u;
    int    n;

    for (n=0; n<N; n++) {
        // Generate "random" variates
        u = (n+0.5) / N;

        // extra n*1e-100 is to prevent compiler
        // optimising for fixed lam
        x = poissinv_v(u, lam+n*1e-100);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Single precision version
void normcdfinv_test_float_scalar(int M) {
    float x, u;
    int m;

    for (m=0; m<M; m++) {
        // Generate "random" variates
        u = (m+0.5) / M;

        x = norminvf(u);

        // needed to prevent compiler discarding everything
        if (x==-999.0) printf("negative x\n");
    }
}

// Single precision version CPU/MIMD-style version
void poissinv_test_float_scalar(int N, float lam) {
    float x, u;
    int n;

    for (n=0; n<N; n++) {
        // Generate "random" variates
        u = (n+0.5f) / N;

        // extra n*1e-20 is to prevent compiler
        // optimising for fixed p
        x = poissinvf(u, lam + n*1e-20f);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

// Single precision version GPU/SIMD-style version
void poissinv_test_float_vector(int N, float lam) {
    float x, u;
    int n;

    for (n=0; n<N; n++) {
        // Generate "random" variates
        u = (n+0.5f) / N;

        // extra n*1e-20 is to prevent compiler
        // optimising for fixed p
        x = poissinvf_v(u, lam + n*1e-20f);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

//
// main code
//

int main(int argc, char **argv) {
    double timer, elapsed;  // timer variable and elapsed time
    float lam;
    int    N, pass, count;

    // number of loop iterations in timing tests
#ifndef COUNT_LAMBDA
#define COUNT_LAMBDA 16
#endif

    char filename[128];
    FILE *fp;

    if (argc > 2) {
        sprintf(filename, "%s", argv[2]);
    } else {
        sprintf(filename, "poissinv_timing.txt");
    }
    fp = fopen(filename, "w");

#define REPEAT 4

    N = (1<<25);

    // execute code

    for (pass=0; pass<4; pass++) {
        // Option to test only double prec or single prec codes

#ifdef DOUBLE
        if (pass % 2 == 0) continue;
#endif
#ifdef SINGLE
        if (pass % 2 == 1) continue;
#endif
#ifndef VECTOR
        if (pass > 1) break;
#endif
#ifdef SKIP_SCALAR
        if (pass < 2) continue;
#endif

        // default parameter values
        lam = 0.125f;

        // change lambda based on input if desired
        if (argc > 1) {
            lam = strtod(argv[1], NULL);
        }

        if (pass==0) {
            PRINTF2(fp, "\nscalar single precision algorithm performance tests (CPU) \n");
        } else if (pass==1) {
            PRINTF2(fp, "\nscalar double precision algorithm performance tests (CPU) \n");
        } else if (pass==2) {
            PRINTF2(fp, "\nvector single precision algorithm performance tests (CPU) \n");
        } else {
            PRINTF2(fp, "\nvector double precision algorithm performance tests (CPU) \n");
        }
        PRINTF2(fp, "---------------------------------- \n");
        PRINTF2(fp, "    lambda   execution time   samples/sec \n");

        for (count=0; count<COUNT_LAMBDA; count++) {
            elapsed_time(&timer);  // initialise timer

            // average over REPEAT runs
            for (int i=0; i<REPEAT; i++) {
                if (pass == 0) {
                    poissinv_test_float_scalar(N, lam);
                } else if (pass==1) {
                    poissinv_test_double_scalar(N, (double) lam);
                } else if (pass==2) {
                    poissinv_test_float_vector(N, lam);
                } else {
                    poissinv_test_double_vector(N, (double) lam);
                }
            }

            elapsed = elapsed_time(&timer);

            if (count>0) { // skip first one (cache effects?)
                PRINTF2(fp, "   %6g      %9.4f     %10.2e \n",
                        lam, elapsed, N*REPEAT/elapsed);
            }
            // update parameter
            /* lam += .25; */
            lam *= 2.0f;
        }

        /* norminv */
        /* Run once first (cache effects?) */
        if (pass == 1 || pass == 3) {
            normcdfinv_test_double_scalar(N);
        } else {
            normcdfinv_test_float_scalar(N);
        }
        elapsed_time(&timer);   // initialise timer
        for (int i=0; i<REPEAT; i++) {
            if (pass == 1 || pass == 3) {
                normcdfinv_test_double_scalar(N);
            } else {
                normcdfinv_test_float_scalar(N);
            }
        }
        elapsed = elapsed_time(&timer);
        PRINTF2(fp, "   normcdfinv  %9.4f     %10.3g \n",
                elapsed, REPEAT*N/elapsed);
    }

    fclose(fp);
    return 0;
}
