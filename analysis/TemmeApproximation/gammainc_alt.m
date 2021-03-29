function S = gammainc_alt(X,A,flag)
    %
    % Incomplete regularized gamma function, more accurate than native
    % MATLAB for large arguments. Computation using Temme uniform
    % asymptotic expansion, see 
    % Chapter 8 Numerical Methods for Special Functions
    % Amparo Gil , Javier Segura and Nico M. Temme
    % https://doi.org/10.1137/1.9780898717822
    % 
    % gammainc_alt(x,a,0) = gammainc(x,a,'lower')
    % gammainc_alt(x,a,1) = gammainc(x,a,'upper')
    %
    if nargin < 3
        flag = 0;
    end
    
    Ai = 1.0 ./ A;
    % Use rlog1(z) = z - log(1+z) function, a fast rational
    % minimax approximation for -0.6 < z < 1.6. This is a lot
    % more accurate when z is near 0 (here A near X).
    Eta = sqrt(2.0*rlog1((X-A).*Ai));
    ix = (X < A);
    Eta(ix) = -Eta(ix);

    B1 =  8.0995211567045583e-16;               S = B1;
    B0 = -1.9752288294349411e-15;               S = B0 + S.*Eta;
    B1 = -5.1391118342426808e-16 + 25.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  2.8534893807047458e-14 + 24.0*B0.*Ai; S = B0 + S.*Eta;
    B1 = -1.3923887224181616e-13 + 23.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  3.3717632624009806e-13 + 22.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =  1.1004392031956284e-13 + 21.0*B1.*Ai; S = B1 + S.*Eta;
    B0 = -5.0276692801141763e-12 + 20.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =  2.4361948020667402e-11 + 19.0*B1.*Ai; S = B1 + S.*Eta;
    B0 = -5.8307721325504166e-11 + 18.0*B0.*Ai; S = B0 + S.*Eta;
    B1 = -2.5514193994946487e-11 + 17.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  9.1476995822367933e-10 + 16.0*B0.*Ai; S = B0 + S.*Eta;
    B1 = -4.3820360184533521e-09 + 15.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  1.0261809784240299e-08 + 14.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =  6.7078535434015332e-09 + 13.0*B1.*Ai; S = B1 + S.*Eta;
    B0 = -1.7665952736826086e-07 + 12.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =  8.2967113409530833e-07 + 11.0*B1.*Ai; S = B1 + S.*Eta;
    B0 = -1.8540622107151585e-06 + 10.0*B0.*Ai; S = B0 + S.*Eta;
    B1 = -2.1854485106799979e-06 +  9.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  3.9192631785224383e-05 +  8.0*B0.*Ai; S = B0 + S.*Eta;
    B1 = -0.00017875514403292177 +  7.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =  0.00035273368606701921 +  6.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =   0.0011574074074074078 +  5.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =   -0.014814814814814815 +  4.0*B0.*Ai; S = B0 + S.*Eta;
    B1 =    0.083333333333333329 +  3.0*B1.*Ai; S = B1 + S.*Eta;
    B0 =    -0.33333333333333331 +  2.0*B0.*Ai; S = B0 + S.*Eta;
    S  = S ./ (1.0 + B1.*Ai);

    S = S.*exp(-0.5*A.*Eta.*Eta)./sqrt(2*pi*A);
    if flag 
        S = 0.5*erfc(Eta.*sqrt(0.5*A))  + S;
    else
        S = 0.5*erfc(-Eta.*sqrt(0.5*A)) - S;
    end

end


