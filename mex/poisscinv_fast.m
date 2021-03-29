% poisscinv_fast.m Help file for poissinv_fast MEX-file.
%  poissinv_fast.c - Fast evaluation of the Poisson inverse complementary
%  cumulative distribution function.
%
%   X = poisscinv_fast(P,LAMBDA) returns the inverse of the Poisson 
%   complementary cdf with parameter lambda. Since the Poisson distribution
%   is discrete, poisscinv_fast returns the smallest value of X, such that 
%   the Poisson complementary cdf evaluated, at X, equals or exceeds P.
%
%   The size of X is the common size of P and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input. 
%
%   Reference:
%      [1]  M.B. Giles, "Algorithm 955: approximation of the inverse 
%           Poisson cumulative distribution function", ACM Transactions on 
%           Mathematical Software, January 2016, Article No.: 7, 
%           https://doi.org/10.1145/2699466
%
%   https://people.maths.ox.ac.uk/gilesm/codes/poissinv/
%   https://github.com/cbeentjes/poissinv