/*==========================================================
 * poisscinv_fast.c -
 *
 * Computes the inverse Poisson complementary CDF efficiently
 *
 * The calling syntax is:
 *
 *		outMatrix = poisscinv_fast(U, lambda)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"  /*Include MEX library*/

/* Trick to convert math.h library FP_NAN into MEX compatible NaN and
 FP_INFINITE into MEX compatible Inf */
#include <math.h>

#ifdef FP_NAN
#undef FP_NAN
#define FP_NAN mxGetNaN()
#endif

#ifdef FP_INFINITE
#undef FP_INFINITE
#define FP_INFINITE mxGetInf()
#endif

#include "poissinv.h" /*Include Giles Inverse Poisson CDF library*/

/* Reset FP_NAN and FP_INFINITE macros */
#include <math.h>

/* MATLAB provides a preprocessor macro, mwsize, that represents size
 * values for integers, based on the platform. The computational
 * routine declares the size of the array as int. Replace the int
 * declaration for variables n and i with mwsize.*/

/*Computational routine prototypes */
void poisscinv_fast_full(double *U, double *lam, double *Z, mwSize n);
void poisscinv_fast_l(double *U, double *lam, double *Z, mwSize n);
void poisscinv_fast_U(double *U, double *lam, double *Z, mwSize n);

/* Gateway/Main function MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])

{
    const mwSize *dims_lambda, *dims_U; /* Dimensions of the input matrices */
    double *lambda_matrix, *U_matrix;   /* Input matrices */
    const mwSize *dims_out;             /* Dimensions of the output matrix */
    mwSize n_out, n_dim;                /* Number of dimensions and elements in output matrix */
    double *out_matrix;                 /* Output matrix */

    /* check for proper number of arguments */
    if (nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:nrhs","Two input arguments required.\n");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:nlhs","Too many output arguments.\n");
    }

    /* make sure the input arguments are type double */
    if (!mxIsDouble(prhs[0])
        || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:notDouble","Input P must be real of type double.\n");
    }
    if (!mxIsDouble(prhs[1])
        || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:notDouble","Input lambda must be real of type double.\n");
    }

    /* check that dimensions are the same for input arguments or scalar */
    if (mxGetNumberOfDimensions(prhs[0])!=mxGetNumberOfDimensions(prhs[1])
	    && mxGetNumberOfElements(prhs[1])!=1
        && mxGetNumberOfElements(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:InputSizeMismatch","Requires non-scalar arguments to match in size.\n");
    }

    /* get the number of dimensions of the inputs */
    dims_lambda = mxGetDimensions(prhs[1]);
    dims_U      = mxGetDimensions(prhs[0]);

    /* check whether the inputs have the same dimensions */
    if (mxGetNumberOfElements(prhs[0])!=1
        && mxGetNumberOfElements(prhs[1])!=1)
    {
        for (mwSize d=0; d<mxGetNumberOfDimensions(prhs[0]); d++) {
            if (dims_lambda[d]!=dims_U[d]) {
                mexErrMsgIdAndTxt("MyToolbox:poisscinv_fast:InputSizeMismatch","Inputs must have the same dimensions.\n");
            }
        }
    }

    /* create a pointer to the real data in the input rate matrix  */
    lambda_matrix = mxGetPr(prhs[1]);

    /* create a pointer to the real data in the input probability matrix  */
    U_matrix = mxGetPr(prhs[0]);

    /* get dimensions of the output matrix */
    if (mxGetNumberOfElements(prhs[0])!=1) {
        n_dim     = mxGetNumberOfDimensions(prhs[0]);
        dims_out    = mxGetDimensions(prhs[0]);
    } else {
        n_dim  = mxGetNumberOfDimensions(prhs[1]);
        dims_out = mxGetDimensions(prhs[1]);
    }

    /* create the output matrix */
    plhs[0] = mxCreateNumericArray(n_dim,dims_out,mxDOUBLE_CLASS,mxREAL);
    n_out   = mxGetNumberOfElements(plhs[0]);

    /* get a pointer to the real data in the output matrix */
    out_matrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    if (mxGetNumberOfElements(prhs[0])==1) {
       poisscinv_fast_U(U_matrix, lambda_matrix, out_matrix, n_out);
    } else if (mxGetNumberOfElements(prhs[1])==1) {
       poisscinv_fast_l(U_matrix, lambda_matrix, out_matrix, n_out);
    } else {
       poisscinv_fast_full(U_matrix, lambda_matrix, out_matrix, n_out);
    }
}

/* Computational routines */
/* Matrix rate lambda and U */
void poisscinv_fast_full(double *U, double *lam, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = poisscinv(U[i],lam[i]);
    }
}

/* Scalar rate lambda */
void poisscinv_fast_l(double *U, double *lam, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = poisscinv(U[i],lam[0]);
    }
}

/* Scalar U */
void poisscinv_fast_U(double *U, double *lam, double *Z, mwSize n)
{
    mwSize i;
    for (i=0; i<n; i++) {
        Z[i] = poisscinv(U[0],lam[i]);
    }
}