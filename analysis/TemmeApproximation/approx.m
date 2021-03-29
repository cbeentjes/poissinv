%
% this routine produces various polynomial approximations to the Temme
% uniform expansion and verifies the accuracy of said approximations.
% Additionaly it verifies the error bounds (delta) for the Temme uniform
% expansion in various scenarios
%

close all;
clear all;
addpath('../../mex')
if ~isfile('../../mex/poissinv_fast.mexa64')
    mex ../../mex/poissinv_fast.c -outdir ../../mex -I../../src/Serial
end

%
% set extent of central region in which approximations are defined
%

rhi = 3.25;
rlo = 0.4;

slo = sqrt(2*(1-rlo + rlo.*log(rlo))) .* sign(rlo-1);
shi = sqrt(2*(1-rhi + rhi.*log(rhi))) .* sign(rhi-1);

file = ' ';
file = char(file,'//  use polynomial approximations in central region',' ');
file = char(file,sprintf('    if ( (s>%10.7gf) && (s<%8.7gf) ) {',slo,shi));
disp(file)

%
% plot f(r)
%

r = linspace(0.0001, rhi);
y = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
plot(r,y)
xlabel('r'); ylabel('f(r)')

% print('-deps2','approx1.eps')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% first, approximate f^{-1} - 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_chebyshev = 200; % number of chebyshev nodes for r

th  = pi*(0.5 + 0:n_chebyshev-1)'/n_chebyshev;
r   = 0.5*(rhi+rlo) - 0.5*(rhi-rlo)*cos(th);
y   = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
s   = (r-1)./y;

%
% least-squares polynomial fit with diagonal weighting
%

degree = 14;

D = diag(1./s);
A = zeros(length(y),degree+1);
for deg = 0:degree
  A(:,degree+1-deg) = y.^deg;
end

A = A(:,2:end-2);
p1 = ((D*A)\(D*(s-1-y/6)))'; % subtract out leading order behaviour
p1 = [p1 1/6 1];             % add back in leading order behaviour

%
% output polynomial approximation
%

file = ' ';
file = char(file,'//  polynomial approximation to f^{-1}(s) - 1',' ');
file = char(file,sprintf('      rm =  %16.9gf;',p1(1)));

for deg = 2:degree-2
  file = char(file,sprintf('      rm =  %16.9gf + rm*s;',p1(deg)));
end
file = char(file,sprintf('      Rm = %16.16g + rm*s;',p1(deg+1)));
file = char(file,'      Rm =              S + S*(Rm*S);');

disp(file)


% Evaluate most terms in single precision is fine
for deg = 1:length(p1)
    p1(deg) = str2double(sprintf('%16.9g',p1(deg)));
end
p1 = double(single(p1));
% Leading order should be done in double precision
p1(end-1) = 1/6;

s2 = polyval(p1,y);
p1 = [p1 0];

%
% plot error
%

figure(2)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

subplot(1,2,1)
plot(y,y.*(s-s2))
xlabel('s'); ylabel('error'); title('    f^{-1} - 1 approximation')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% second, approximate c_0(r)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_chebyshev = 200; % number of chebyshev nodes for r

th  = pi*(0.5 + 0:n_chebyshev-1)'/n_chebyshev;
r   = 0.5*(rhi+rlo) - 0.5*(rhi-rlo)*cos(th);
y   = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
c   = log(y.*sqrt(r)./(r-1)) ./ log(r);

%
% least-squares polynomial fit with diagonal weighting
%

degree = 12;

D = diag(1./r.^0.5);
A = zeros(length(r),degree);
for deg = 1:degree
  A(:,degree+1-deg) = (r-1).^deg;
end

p2 = ((D*A)\(D*(c-1/3)))'; % subtract out leading order behaviour
p2 = [p2 1/3];             % add back in leading order behaviour

%
% output polynomial approximation
%

file = ' ';
file = char(file,'//  polynomial approximation to correction c0(r)',' ');
file = char(file,sprintf('      t  = %16.9gf;',p2(1)));

for deg = 1:degree
  file = char(file,sprintf('      t  = %16.9gf + t*rm;',p2(deg+1)));
end

disp(file)

% Evaluate the polynomial in single precision
for deg = 1:length(p2)
    p2(deg) = str2double(sprintf('%16.9g',p2(deg)));
end
p2 = double(single(p2));

c2 = polyval(p2,r-1);

%
% plot error
%
subplot(1,2,2)
plot(r,c-c2); 
xlabel('r'); ylabel('error'); title('    c_0 approximation')

% print('-deps2','approx2.eps')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% third, additional O(1/lam) correction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xs = [10 10.5 11 12 13 16 20 50 80 100];

% Increase domain of approximation slightly to incorporate the effect that
% lam_min, lam_max selection below is marginally incorrect at small values
% of x.
rhi = 1.009*rhi;

degree = 10;

D = eye(1000*length(xs));
A = zeros(1000*length(xs),degree+1);
c = zeros(1000*length(xs),1);

for k = 1:length(xs)
  x = xs(k);

  lam_min = x/rhi;
  lam_max = x/rlo;  

  n_nodes = 1e3;
  lam = linspace(lam_min,lam_max,n_nodes);
  w   = zeros(size(lam));
  x3  = zeros(size(lam));
  
%
% for accuracy in tails, need to select correct variant
%

  upper = find(x<lam);
  w(upper) = norminv(gammainc(lam(upper),x,'upper'));
  lower = find(x>=lam);
  w(lower) = -norminv(gammainc(lam(lower),x,'lower'));

%
% use Newton iteration to find rm
%

  s  = w./sqrt(lam);
  rm = finv(s);
  c0 = log(s.*sqrt(1+rm)./rm) ./ log(1+rm);
  x2 = lam + lam.*rm;

%
% full correction of leading order error
%
  t = abs(rm)<0.0001;
  x3(t)  = x2(t) + 1/3 - rm(t)/36;
  x3(~t) = x2(~t) + c0(~t);

  c(n_nodes*(k-1)+1:n_nodes*k) = (x - x3)';
  D(n_nodes*(k-1)+1:n_nodes*k,n_nodes*(k-1)+1:n_nodes*k) = ...
       diag((1:n_nodes).^(-0.75).*(n_nodes:-1:1).^(-0.25));

  for deg = 0:degree
    A(n_nodes*(k-1)+1:n_nodes*k,degree+1-deg) = rm.^deg ./ lam;
  end
end

%
% compute polynomial approximation
%

p3 = ((D*A)\(D*c))';

file = ' ';
file = char(file,'//  O(1/lam) correction',' ');
file = char(file,sprintf('      x  = %16.9gf;',p3(1)));

for deg = 1:degree
  file = char(file,sprintf('      x  = %16.9gf + x*rm;',p3(deg+1)));
end

disp(file)

% Return to usual domain of approximation
rhi = rhi/1.009;

% Evaluate the polynomial in single precision
for deg = 1:length(p3)
    p3(deg) = str2double(sprintf('%16.9g',p3(deg)));
end
p3 = double(single(p3));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 1a
%
% Error for a fixed value of x, while varying lambda, for the polynomial
% approximation to Temme uniform expansion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

xs  = [10 3.3e4 exp(linspace(log(1e1),log(1e9),1e3))];
err = zeros(size(xs));
err_s = zeros(size(xs));

for k = 1:length(xs)
  x = xs(k);
  
% single precision range
  lam_min_s = max([gammaincinv(1e-7, x,'lower') x/rhi 4]);
  lam_max_s = min([gammaincinv(1e-38,x,'upper') x/rlo]);  
  
% double precision range
  lam_min = max([gammaincinv(1e-16, x,'lower') x/rhi 4]);
  lam_max = min([gammaincinv(1e-310,x,'upper') x/rlo]);
 
  % Account for the fact that for small x the bounds x/rhi & x/rlo are not
  % quite appropriate. lam_min and lam_max are too large due to higher
  % order corrections.
  if lam_min == x/rhi
      fac     = 2;
      for kount = 1:20
         while(gammainc(lam_min/fac,x,'lower')>0 && -norminv(gammainc(lam_min/fac,x,'lower'))/sqrt(lam_min/fac) < shi)
            lam_min = lam_min/fac;
         end
         fac = sqrt(fac);
      end
      lam_min = max(4,lam_min);  
  end    
  if lam_max == x/rlo
      fac     = 1.2;      
      for kount = 1:20
         while(gammainc(lam_max/fac,x,'upper')>0 && norminv(gammainc(lam_max/fac,x,'upper'))/sqrt(lam_max/fac) < slo)
            lam_max = lam_max/fac;
         end
         fac = sqrt(fac);
      end          
  end    
  

  lam = linspace(lam_min,lam_max,1e3);
  
  % Index for single precision range
  ix_s = and(lam >= lam_min_s, lam <= lam_max_s);
  
  w   = zeros(size(lam));

  lower = find(x>=lam);
  upper = find(x< lam);  
  
  % Built-in MATLAB gammainc function becomes too inaccurate for large
  % argument values (relative error ~1e-9 rather than machine precision)
  if lam_min > 1e3
      w(lower) = -norminv(gammainc_alt(lam(lower),x,0));      
      w(upper) =  norminv(gammainc_alt(lam(upper),x,1));
  else
      w(lower) = -norminv(gammainc(lam(lower),x,'lower'));
      w(upper) =  norminv(gammainc(lam(upper),x,'upper'));      
  end
  
  s = w./sqrt(lam);
  
  rm = polyval(p1,s);

  x2 = lam.*(1+rm) + polyval(p2,rm) + polyval(p3,rm)./lam;

  if k<=2
    subplot(1,2,k)
    plot(lam,x2-x)
    title(['$x = ' num2str(x) '$'],'interpreter','latex')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error')
  else
    err(k) = max(abs(x2-x));
    err_s(k) = max(abs(x2(ix_s)-x));
  end
end

% print('-deps2','approx3.eps')

figure(4)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
semilogx(xs(3:end),err(3:end));
hold on
semilogx(xs(3:end),err_s(3:end),'r--');
legend('FP64','FP32')
xlabel('$x$','interpreter','latex')
ylabel('maximum error')
%axis([1e1 4e3 0 1e-4])

% print('-deps2','approx4.eps')

% max(err)
% max(err_s)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 1b
%
% Error for a fixed value of lambda, while varying x, for the polynomial
% approximation to Temme uniform expansion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag to determine whether to focus on CPU or GPU version algorithm
cpu_flag = false;
    
figure(5)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

lambda  = [10 1e5 exp(linspace(log(4),log(1e9),1e2))];
err   = zeros(numel(lambda),2);
err_s = zeros(numel(lambda),2);

min_x = zeros(numel(lambda),1);
min_x_s = zeros(numel(lambda),1);

for k = 1:length(lambda)
  lam = lambda(k);
  
% double precision range
  x_min = poissinv_fast(1e-310,lam);
  x_max = poissinv_fast(1-1e-16,lam);

% single precision range
  x_min_s = poissinv_fast(1e-38,lam);
  x_max_s = poissinv_fast(1-1e-7,lam);
  
  x_min   = max(x_min,10);
  x_min_s = max(x_min_s,10);

  X = linspace(x_min,x_max,1000);
  % Index for single precision range
  ix_s = and(X >= x_min_s, X <= x_max_s);    
  
  w = zeros(size(X));
  
  upper = find(X<lam);
  lower = find(X>=lam);
  
  % Built-in MATLAB gammainc function becomes too inaccurate for large
  % argument values (relative error ~1e-9 rather than machine precision)
  if lam > 1e3
      w(lower) = -norminv(gammainc_alt(lam,X(lower),0));
      w(upper) =  norminv(gammainc_alt(lam,X(upper),1));
  else
      w(lower) = -norminv(gammainc(lam,X(lower),'lower'));      
      w(upper) =  norminv(gammainc(lam,X(upper),'upper'));      
  end
  
  % Remove central CPU region
  if cpu_flag
    X = X(abs(w)>3);
    w = w(abs(w)>3);
  end
  
  s = w./sqrt(lam);
  
  valid = and(s > slo, s < shi);
  w = w(valid);
  X = X(valid);
  s = s(valid);
  ix_s = ix_s(valid);
  
  rm = polyval(p1,s);

  x2 = lam.*(1+rm) + polyval(p2,rm) + polyval(p3,rm)./lam;
  
  if cpu_flag
      % Asymptotic expansion for when s close to zero
      rm = 0;
      rm =  19/34020 + rm.*s;
      rm = -23/17280 + rm.*s;
      rm =  1/270 + rm.*s;
      rm = -1/72  + rm.*s;
      rm =  1/6   + rm.*s;
      rm = 1      + rm.*s;
      rm = rm.*s;

      c0 = 0;
      c0 = -137/38880 + c0.*s;
      c0 = 7/810 + c0.*s;
      c0 = -1/36 + c0.*s;
      c0 = 1/3 + c0.*s;

      x3 = lam + (lam*rm + c0);

      % Asymptotic expansion in single precision when s close to zero
      rm = 0;
      rm =  1/6  + rm.*s;
      rm = 1     + rm.*s;
      rm = rm.*s;  

      c0 = 0;    
      c0 = 1/3 + c0.*s;

      x3_s = lam + (lam*rm + c0);  
  end
  
  min_x(k)   = x_min;
  min_x_s(k) = x_min_s;

  if k<=2
    subplot(1,2,k)
    plot(X,x2-X)
    title(['$\lambda = ' num2str(lam) '$'],'interpreter','latex')
    xlabel('$x$','interpreter','latex'); ylabel('error')
  else
    err(k,1)   = max(abs(x2-X));
    err_s(k,1) = max(abs(x2(ix_s)-X(ix_s)));
    
    if cpu_flag
        idx = abs(s) < 1e-2;
        if any(idx)
            err(k,2) = max(abs(x3(idx)-X(idx)));
            err_s(k,2) = max(abs(x3_s(idx)-X(idx)));
        end
    end
  end
end

% print('-deps2','approx3b.eps')

figure(6)
loglog(lambda(3:end),err(3:end,:));
hold on 
loglog(lambda(3:end),err_s(3:end,:));
loglog(lambda(3:end),eps(single(min_x_s(3:end))))

loglog(lambda(3:end),2.5e-2./lambda(3:end),'--')
legend('FP64 (Newton)','FP64 (Expansion)','FP32 (Newton)','FP32 (Expansion)',...
    'Single precision relative precision','0.025/\lambda','location','best')
ylabel('maximum error')
xlabel('$\lambda$','interpreter','latex')
%axis([1e1 4e3 0 1e-4])

% print('-deps2','approx4b.eps')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 2a
%
% Error for a fixed value of x, while varying lambda, for the Newton
% iteration to solve Temme uniform expansion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag to determine whether to focus on CPU or GPU version algorithm
cpu_flag = true;

figure(7)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

xs  = [1e1 1e4 exp(linspace(log(10),log(1e8),1000))];
err1 = zeros(size(xs));
err2 = zeros(size(xs));
err1_s = zeros(size(xs));
err2_s = zeros(size(xs));

for k = 1:length(xs)
  x  = xs(k);

% double precision range

% brute force search for lam_min because gammaincinv fails

  lam_min = x;
  fac     = 2;
  for kount = 1:20
     while(gammainc(lam_min/fac,x,'lower')>1e-310)
       lam_min = lam_min/fac;
     end
     fac = sqrt(fac);
  end
  lam_min = max(4,lam_min);

  lam_max = gammaincinv(1e-310,x,'upper');

% single precision range
  lam_min_s = max(4, gammaincinv(6.e-8,x,'lower'));
  lam_max_s =        gammaincinv(1e-38,x,'upper');

  lam = linspace(lam_min,lam_max,1000);
  
  % Index for single precision range
  ix_s = and(lam >= lam_min_s, lam <= lam_max_s);
  
  w   = zeros(size(lam));

  lower = find(x>=lam);
  upper = find(x< lam);  
  
  % Built-in MATLAB gammainc function becomes too inaccurate for large
  % argument values (relative error ~1e-9 rather than machine precision)
  if lam_min > 1e3
      w(lower) = -norminv(gammainc_alt(lam(lower),x,0));      
      w(upper) =  norminv(gammainc_alt(lam(upper),x,1));
  else
      w(lower) = -norminv(gammainc(lam(lower),x,'lower'));
      w(upper) =  norminv(gammainc(lam(upper),x,'upper'));      
  end  
  
  % Remove central CPU region (done by Normal approximation)
  if cpu_flag
    lam = lam(abs(w)>3);
    w = w(abs(w)>3);
    ix_s = ix_s(abs(w)>3);
  end  

  s  = w./sqrt(lam);
  
  % Remove asymptotic expansion when s close to zero when CPU flag active
  if cpu_flag
    lam = lam(abs(s)>1e-2);
    w   = w(abs(s)>1e-2);
    s   = s(abs(s)>1e-2);
    ix_s= ix_s(abs(s)>1e-2);
   % Remove polynomial approximation for slo < s < shi when GPU flag active
  else
    valid = or(s <= slo, s >= shi);
    w = w(valid);
    lam = lam(valid);
    s = s(valid);
    ix_s = ix_s(valid);
  end
  
  rm = finv(s);
  c0 = log(s.*sqrt(1+rm)./rm) ./ log(1+rm);
  x1 = lam+lam.*rm;
  
  x1 = x1 + c0;

  x2 = x1 - 0.0218./(x1+0.065*lam);

  if k<=2
    subplot(2,2,k)
    plot(lam,x1-x,'.')
    title(['x = ' num2str(x)])
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 1')
    subplot(2,2,k+2)
    plot(lam,x2-x,'.')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 2')
  else
      if ~isempty(x1)
        err1(k) = max(abs(x1-x));
        err2(k) = max(abs(x2-x));      
      end
      if any(ix_s)
        err1_s(k) = max(abs(x1(ix_s)-x));
        err2_s(k) = max(abs(x2(ix_s)-x)); 
      end
  end
end

% print('-deps2','approx5.eps')

figure(8)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
loglog(xs(3:end),err1(3:end),'-k',xs(3:end),err2(3:end),'-.k');
hold on
loglog(xs(3:end),err1_s(3:end),'-b',xs(3:end),err2_s(3:end),'-.b');

loglog(xs(3:end),0.01./xs(3:end),'r')
loglog(xs(3:end),0.002./xs(3:end),'r--')

xlabel('$x$','interpreter','latex')
ylabel('maximum error')
legend('error 1 (FP64)','error 2 (FP64)','error 1 (FP32)','error 2 (FP32)','0.01/x','0.002/x')

% print('-deps2','approx6.eps')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 2b
%
% Error for a fixed value of lambda, while varying x, for the Newton
% iteration to solve Temme uniform expansion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag to determine whether to focus on CPU or GPU version algorithm
cpu_flag = true;

figure(9)
% pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

lambda  = [10 7e2 exp(linspace(log(1e1),log(1e8),1e2))];
err1 = zeros(numel(lambda),1);
err2 = zeros(numel(lambda),1);
err1_s = zeros(numel(lambda),1);
err2_s = zeros(numel(lambda),1);

for k = 1:length(lambda)
  lam  = lambda(k);

% double precision range

  x_min = poissinv_fast(1e-310,lam);
  x_max = poissinv_fast(1-1e-16,lam);

% single precision range
  x_min_s = poissinv_fast(1e-38,lam);
  x_max_s = poissinv_fast(1-1e-7,lam);

  x_min = max(10, x_min);
  x_min_s = max(10, x_min_s);

  X = linspace(x_min,x_max,1000);

  % Index for single precision range
  ix_s = and(X >= x_min_s, X <= x_max_s);    
  
  w = zeros(size(X));
  
  upper = find(X<lam);
  lower = find(X>=lam);
  
  % Built-in MATLAB gammainc function becomes too inaccurate for large
  % argument values (relative error ~1e-9 rather than machine precision)
  if lam > 1e5
      w(lower) = -norminv(gammainc_alt(lam,X(lower),0));
      w(upper) =  norminv(gammainc_alt(lam,X(upper),1));
  else
      w(lower) = -norminv(gammainc(lam,X(lower),'lower'));      
      w(upper) =  norminv(gammainc(lam,X(upper),'upper'));      
  end  

  % Remove central CPU region (done by Normal approximation)
  if cpu_flag
    X = X(abs(w)>3);
    w = w(abs(w)>3);
  end   
  
  s = w./sqrt(lam);
    
  % Remove asymptotic expansion when s close to zero when CPU flag active
  if cpu_flag
    X   = X(abs(s)>1e-2);
    w   = w(abs(s)>1e-2);
    s   = s(abs(s)>1e-2);
    ix_s= ix_s(abs(s)>1e-2);
   % Remove polynomial approximation for slo < s < shi when GPU flag active
  else
    valid = or(s <= slo, s >= shi);
    w = w(valid);
    X = X(valid);
    s = s(valid);
    ix_s = ix_s(valid);
  end
  
  
  rm = finv(s);
  c0 = log(s.*sqrt(1+rm)./rm) ./ log(1+rm);
  x1 = lam.*rm;

  x1 = x1  + c0;

  x2 = x1 - 0.0218./(x1+1.065*lam);
  
  x1 = x1 + lam;
  x2 = x2 + lam;
  
  if k<=2
    subplot(2,2,k)
    plot(X,x1-X)
    title(['\lambda = ' num2str(lam)])
    xlabel('$x$','interpreter','latex'); ylabel('error 1')
    subplot(2,2,k+2)
    plot(X,x2-X)
    xlabel('$x$','interpreter','latex'); ylabel('error 2')
  else
      if ~isempty(X)
        err1(k) = max(abs(x1-X));
        err2(k) = max(abs(x2-X));
      end
      if ~isempty(X(ix_s))
        err1_s(k) = max(abs(x1(ix_s)-X(ix_s)));
        err2_s(k) = max(abs(x2(ix_s)-X(ix_s)));
      end
  end
end


% print('-deps2','approx5b.eps')

figure(10)
% pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
loglog(lambda(3:end),err1(3:end),'-k',lambda(3:end),err2(3:end),'-.k');
hold on
loglog(lambda(3:end),err1_s(3:end),'-b',lambda(3:end),err2_s(3:end),'-.b');
xlabel('$\lambda$','interpreter','latex')
ylabel('maximum error')
legend('error 1 (FP64)','error 2 (FP64)','error 1 (FP32)','error 2 (FP32)')


% print('-deps2','approx6b.eps')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Newton iteration to calculate f^{-1}(s) - 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rm = finv(s)

    t     = (s > -sqrt(2));
    r(~t) = 0;
    r(t)  = max(0.1,1+s(t));
    t     = t & (s~=0);

    while (sum(t)>0)
      r_old = r;
      logr  = log(r(t));
      f = sqrt(2*(1 - r(t) + r(t).*logr)) .* sign(r(t)-1);
      r(t) = r_old(t) - (f-s(t)).*f./logr;
      r(t) = max(r(t), 0.1*r_old(t));
      t = abs(r-r_old)>1e-5;
    end

    rm = r-1;

end
