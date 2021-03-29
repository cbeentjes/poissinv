//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>

//
// CUDA header file
//

#include <cuda.h>

//
// my header file
//

#include "poissinv_cuda.h"

// number of samples per thread
#define SAMPLES 100.0

__global__ void normcdfinvf_test(int M) {

    float x, u;
    int   tid = threadIdx.x + blockIdx.x*blockDim.x;

    u = (tid + 0.5f) / M;

    for (int n=0; n<SAMPLES; n++) {
        u += 1e-10f;

        x = normcdfinvf(u);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

__global__ void poissinvf_test(int M, float lam) {

    float x, u;
    int   tid = threadIdx.x + blockIdx.x*blockDim.x;

    u = (tid + 0.5f) / M;

    for (int n=0; n<SAMPLES; n++) {
        u += 1e-10f;

        // Mixed case used when lam < 0
        // lam takes values from {1,2,4,8,16,32,64,128}
        if (lam<0.0f) {
            int n = 1 << (tid & 7);
            lam   = 1.0f * (float) n;
        }
        x = poissinvf(u, lam);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

__global__ void normcdfinv_test(int M) {

    float x, u;
    int   tid = threadIdx.x + blockIdx.x*blockDim.x;

    u = (tid + 0.5f) / M;

    for (int n=0; n<SAMPLES; n++) {
        u += 1e-10f;

        x = normcdfinv((double) u);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}

__global__ void poissinv_test(int M, float lam) {

    float x, u;
    int   tid = threadIdx.x + blockIdx.x*blockDim.x;

    u = (tid + 0.5f) / M;

    for (int n=0; n<SAMPLES; n++) {
        u += 1e-10f;

        // Mixed case used when lam < 0
        // lam takes values from {1,2,4,8,16,32,64,128}
        if (lam<0.0f) {
            int n = 1 << (tid & 7);
            lam   = 1.0f * (float) n;
        }
        x = poissinv((double) u, (double) lam);

        // needed to prevent compiler discarding everything
        if (x==-999.0f) printf("negative x\n");
    }
}



//
// main code
//

int main(int argc, char **argv) {
    float lam;
    int   M, nblocks, nthreads;

    // number of loop iterations in timing tests
#ifndef COUNT_LAMBDA
#define COUNT_LAMBDA 16
#endif

    // CUDA timing

    float milli;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // set number of blocks, and threads per block

    M = (1<<27);
    nthreads = 256;
    nblocks  = M / nthreads;

    // execute kernels

    for (int pass=0; pass<2; pass++) {
        // Option to test only double prec or single prec codes
#ifdef DOUBLE
        if (pass % 2 == 0) continue;
#endif
#ifdef SINGLE
        if (pass % 2 == 1) continue;
#endif

        if (pass==0)
            printf("\nsingle precision performance tests (GPU)\n");
        else
            printf("\ndouble precision performance tests (GPU)\n");
        printf("---------------------------------- \n");
        printf("  lambda   execution time   samples/sec \n");

        // default parameter values
        lam = 0.125f;

        // change lambda based on input if desired
        if (argc > 1) {
            lam = strtod(argv[1], NULL);
        }

        /* Fixed parameters */
        for (int count=0; count<=COUNT_LAMBDA; count++) {

            cudaEventRecord(start);

            if (pass==0)
                poissinvf_test<<<nblocks,nthreads>>>(M, lam);
            else
                poissinv_test<<<nblocks,nthreads>>>(M, lam);

            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&milli, start, stop);

            // factor SAMPLES due to repeat in test routines
            // factor 1e3 due to timing in milliseconds
            if (count>0) { // skip first one for more accurate timing (cache effects?)
                printf("   %6g      %9.4f     %10.3g \n",
                        lam, milli, float(M)*SAMPLES*1e3/milli);
            }
#ifdef PLUS
            lam += 1.0f;
#else
            lam *= 2.0f;
#endif
        }

        /* Mixed parameters */
        // Run once first for more accurate timing
        if (pass == 0)
            poissinvf_test<<<nblocks,nthreads>>>(M, -1.0f);
        else
            poissinv_test<<<nblocks,nthreads>>>(M, -1.0f);

        cudaEventRecord(start);
            if (pass == 0)
                poissinvf_test<<<nblocks,nthreads>>>(M, -1.0f);
            else
                poissinv_test<<<nblocks,nthreads>>>(M, -1.0f);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milli, start, stop);
        printf("    mixed      %9.4f     %10.3g \n",
                milli, float(M)*SAMPLES*1e3/milli);

        /* normcdfinv */
        if (pass == 0)
            normcdfinvf_test<<<nblocks,nthreads>>>(M);
        else
            normcdfinv_test<<<nblocks,nthreads>>>(M);

        cudaEventRecord(start);
        if (pass == 0)
            normcdfinvf_test<<<nblocks,nthreads>>>(M);
        else
            normcdfinv_test<<<nblocks,nthreads>>>(M);

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milli, start, stop);

        printf("\n normcdfinv    %9.4f     %10.3g \n",
                milli, float(M)*SAMPLES*1e3/milli);

    }

    // CUDA exit -- needed to flush printf write buffer

    cudaDeviceReset();
    return 0;
}
