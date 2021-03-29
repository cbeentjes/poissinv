function y = rlog1(x)
%
%   rlog1(x) = x - log(1+x)
%
%   Using minimax approximation for -0.618 < x < 1.618
%   Coefficients from rlog1Expansion.nb notebook.

    % Rational minimax coefficients 
    P0 =   0.2000000000000000;
    P1 =  -0.3636535967811319;
    P2 =   0.2015244511825799;
    P3 =  -0.03274937605228191;
    P4 =   0.00004542775258423288;
    Q1 =  -2.532553698191349;
    Q2 =   2.261033627633934;     
    Q3 =  -0.8263420776370621;
    Q4 =   0.1008870710149048;
    
    r = x./(x+2);
    t = r.*r;

    w = 1/3 + t.*((((P4*t + P3).*t + P2).*t + P1).*t + P0) ./ ...
              ((((Q4*t + Q3).*t + Q2).*t + Q1).*t + 1);
          
    y = t.*(x + 2 - 2*r.*w);           
    
end