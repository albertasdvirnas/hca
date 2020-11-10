function [p] = trunc_normal_ev_pdf(x, m, n, s, a, b)
if s < 0; p = nan; return; end
if nargin < 5
    a = 0;
end
if nargin < 6   
    b = 1;
end

if length(m) > 1 || length(s) > 1
  xs = size(x);
  x = x(:);
  m = m(:)';
  s = s(:)';
end

eta = (x - m)./s;
alpha = (a - m)./s;
beta = (b - m)./s;

npdf = @(x) exp(-.5*x.^2)/sqrt(2*pi);
ncdf = @(x) .5*(1+erf(x/sqrt(2)));

if length(m) > 1 || length(s) > 1
  p = mean((x >= a & x <= b).*(n.*((ncdf(eta) - ncdf(alpha))./(ncdf(beta) - ncdf(alpha))).^(n-1)).*npdf(eta)./(ncdf(beta) - ncdf(alpha))./s, 2);
  p = reshape(p, xs);
else
  p = (x >= a & x <= b).*(n.*((ncdf(eta) - ncdf(alpha))./(ncdf(beta) - ncdf(alpha))).^(n-1)).*npdf(eta)./(ncdf(beta) - ncdf(alpha))./s;
end
p = max(p, 10^-digits);
