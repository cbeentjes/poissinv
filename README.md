# poissinv: Fast evaluation of the inverse Poisson CDF

## Synopsis

The function [normcdfinv(u)](https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html#group__CUDA__MATH__DOUBLE_1g78e93df6c3fbade8628d33e11fc94595) in NVIDIA's [CUDA maths library](https://docs.nvidia.com/cuda/cuda-math-api/index.html) (and similar functions in various libraries for CPU execution) computes the inverse of the cumulative distribution function (CDF) for the Normal distribution. It can be used to convert uniformly distributed pseudo-random (or quasi-random) numbers into psuedo-random (or quasi-random) Normals.

Here we present an update to the software accompanying [TOMS Algorithm 955](https://doi.org/10.1145/2699466). This implements algorithms for a function poissinv(u) which performs the corresponding inversion task for Poisson distributions, and can be used to generate Poisson random variates.

This is mostly based on previous work by M. B. Giles on the [inverse Poission CDF function](https://people.maths.ox.ac.uk/gilesm/codes/poissinv/).

For mathematical details, see the [TOMS paper](https://doi.org/10.1145/2699466)

	@article{giles2016,
	author = {Giles, Michael B.},
	title = {Algorithm 955: Approximation of the Inverse Poisson Cumulative Distribution Function},
	year = {2016},
	volume = {42},
	number = {1},
	doi = {10.1145/2699466},
	journal = {ACM Trans. Math. Softw.},
	month = jan,
	articleno = {7},
	}

For more information on fast inverse CDF computation [click here](https://people.maths.ox.ac.uk/~gilesm/codes/).

#### Main changes relative to original TOMS Algorithm 955

1. Inclusion of a single-precision poissinvf function for CPU use.
1. Improved accuracy; the L<sub>1</sub>-error (see Section 7 of the original TOMS paper) for both single precision and double precision routines grows like &radic;&lambda;, rather than linear in &lambda;, without a performance penalty.

The latter is mainly thanks to:
* Accurate evaluation of the Poisson CDF C(x) so that the error does not grow with &lambda; or x (c.f. Section 4 of the original TOMS paper). This relies on the accurate implementation of rlog1(x) = x - log(1+x) for small x.
* Error bounds &delta; of the fast approximations have been updated. Notably, the Temme approximation used for GPU-computations, Q<sub>T3</sub>, must use the error bound &delta; = 3.3 &#xd7; 10<sup>-6</sup> in double precision.

#### Performance (CPU)

Samples/sec on a single core from a 6-core 35W 2.2Ghz Intel Core i5-9500T, using the compiler flag -O3. 

For comparison see also the number of samples/sec for the inverse normal CDF function. The cost of poissinv(f) is approximately 2-4 times the cost of normcdfinv(f).

| <br>&lambda; | GCC 10.2<br>poissinvf | GCC 10.2<br>poissinv | CLANG 11.0<br>poissinvf | CLANG 11.0<br>poissinv | ICC 19.0<br>poissinvf | ICC 19.0<br>poissinv |
| :--------: | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: |
| 2   | 8.16e07 | 7.02e07 | 8.57e07 | 7.28e07 | 7.86e07 | 6.82e07 |
| 8   | 4.66e07 | 4.33e07 | 4.76e07 | 4.36e07 | 4.54e07 | 4.15e07 |
| 32  | 4.71e07 | 3.43e07 | 4.72e07 | 3.39e07 | 4.40e07 | 3.42e07 |
|128  | 4.75e07 | 3.46e07 | 4.76e07 | 3.43e07 | 4.44e07 | 3.48e07 |
| normcdfinvf | 1.62e08 | | 1.53e08 |  | 1.69e08 | |
| normcdfinv  |  | 9.61e07 |  | 9.29e07 |  | 9.52e07|


#### Performance (GPU/CUDA)

Samples/sec using CUDA 11.2 (code compiled with -O3 and --use_fast_math flags) on three NVIDIA GPUs: 

1. a 250W GeForce RTX 2080 Ti which is a consumer graphics card from 2018 based on the Turing GPU architecture, 
2. a 235W Tesla K40m which is a HPC card from 2013 based on the Kepler architecture, and
3. a 235W Quadro GP100 which is a high-end workstation card from 2016 based on the Pascal architecture. 

For comparison see also the number of samples/sec for the inverse normal CDF function. The cost of poissinv(f) is approximately 2-6 times the cost of normcdfinv(f).

| <br>&lambda; | GeForce RTX 2080 Ti<br>poissinvf | GeForce RTX 2080 <br>poissinv | K40m<br>poissinvf | K40m<br>poissinv | Quadro GP100<br>poissinvf | Quadro GP100<br>poissinv |
| :--------: | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: | :-----------: |
| 2   | 1.99e11 | 4.71e09 | 2.42e10 | 8.26e09 | 6.30e10 | 3.04e10 |
| 8   | 8.16e10 | 1.46e09 | 1.05e10 | 2.88e09 | 3.02e10 | 1.13e10 |
| 32  | 1.14e11 | 2.74e09 | 1.51e10 | 5.21e09 | 4.58e10 | 2.04e10 |
|128  | 1.14e11 | 2.73e09 | 1.51e10 | 5.22e09 | 4.58e10 | 2.03e10 |
|mixed| 7.21e10 | 1.24e09 | 9.01e09 | 2.50e09 | 2.53e10 | 9.76e09 |
| normcdfinvf | 4.22e11 | | 5.43e10 |  | 1.60e11 | |
| normcdfinv  |  | 3.80e09 |  | 1.02e10 |  | 3.59e10|


## Authors

* M. B. Giles   <mike.giles@maths.ox.ac.uk>
* C. H. L. Beentjes <casper.beentjes@maths.ox.ac.uk>

## Installation

To download:

    $ git clone https://github.com/cbeentjes/poissinv.git
    $ cd poissinv

## Matlab integration via MEX files

The fast (serial) C routines can be directly used in MATLAB applications via their MEX interfaces provided. To compile the interfaces run the MATLAB helper file compile_interface.m via the MATLAB IDE
    
    $ cd poissinv/mex
    $ compile_interface

To add inverse CDF and inverse complementary CDF MEX routines to MATLAB's search path via command line:

    $ export MATLABPATH=/path/to/poissinv/mex/:$MATLABPATH

To add said routines via MATLAB IDE:

    $ addpath(/path/to/poissinv/mex/)

To make this permanent consider [adding this line to your startup.m file](mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html).

## Repository content

* [/analysis](https://github.com/cbeentjes/poissinv/tree/main/analysis/) contains MATLAB and Mathematica code used to construct and verify parts of the codebase.

* [/mex](https://github.com/cbeentjes/poissinv/tree/main/mex/) contains MEX wrappers to call the fast inverse CDF and inverse complementary CDF routines as MATLAB routines.

* [/src](https://github.com/cbeentjes/poissinv/tree/main/src/) contains several header files with the core C-routines to evaluate the inverse CDF and complementary CDF functions. 

* [/test](https://github.com/cbeentjes/poissinv/tree/main/test/) contains test code to verify accuracy and test performance of the inverse CDF routines routines. The validation part uses a quad-precision function which requires gcc's [libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) library.
 
## Licensing and acknowledgements

This code is freely available to all under a GPL license -- anyone requiring a more permissive license for commercial purposes should contact M. B. Giles.

We would be grateful if academic users would cite the paper above.
