%%
% Helper file to compile the MEX interfaces for the fast inverse CDF and
% inverse complementary CDF functions for the Poisson distribution

% Double precision variants
mex poissinv_fast.c -I../src/Serial
mex poisscinv_fast.c -I../src/Serial

% Single precision variants
mex poissinvf_fast.c -I../src/Serial
mex poisscinvf_fast.c -I../src/Serial
