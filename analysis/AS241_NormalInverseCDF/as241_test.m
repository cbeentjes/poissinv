%
% This MATLAB routine tests the accuracy of the algorithms in 
%
% Algorithm AS 241: The Percentage Points of the Normal Distribution,
% Wichura, M.J., Applied Statistics, 37(8):477-484, 1988
%

close all
clear all

%
% test main range
%

x1 = linspace(0,0.9,10000);
x1 = x1(2:end);

x1f = single(x1);

y1  = [];
y1f = [];

for p = x1
  y1  = [y1  norminv_as241(p)];
end

for p = x1f
  y1f  = [y1f  norminvf_as241(p)];
end

%
% test tail
%

x2 = 10.^linspace(-310,-1,10000);
x2f = single(10.^linspace(-45,-1,10000));


y2  = [];
y2f = [];

for p = x2
  y2 = [y2 norminv_as241(p)];    
end

for p = x2f
  y2f = [y2f norminvf_as241(p)];
end

%
% plot relative errors
%

figure
subplot(2,1,1)
plot(x1,(y1-norminv(x1))./eps(y1))
title('relative error in main interval (double precision) (ULP)')
xlabel('p')

subplot(2,1,2)
semilogx(x2,(y2-norminv(x2))./eps(y2))
title('relative error in tail (ULP)')
xlabel('p')

figure
subplot(2,1,1)
plot(x1f,(y1f-norminv(x1f))./eps(y1f))
title('relative error in main interval (single precision) (ULP)')
xlabel('p')

subplot(2,1,2)
semilogx(x2f,(y2f-norminv(x2f))./eps(y2f))
title('relative error in tail (ULP)')
xlabel('p')
