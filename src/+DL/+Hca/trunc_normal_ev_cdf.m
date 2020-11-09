function p = trunc_normal_ev_cdf(x, m, n, s, a, b, extraPrecision)
if s < 0; p = nan; return; end
if nargin < 5
  a = -1;
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

ncdf = @(x) .5*(1+erf(x/sqrt(2)));
vpa_ncdf = @(x) .5*(vpa(1)+erf(vpa(x/sqrt(2))));

if nargin < 7 || ~extraPrecision
  if length(m) > 1 || length(s) > 1
    p = mean((x >= a & x <= b).*((ncdf(eta) - ncdf(alpha))./(ncdf(beta) - ncdf(alpha))).^n, 2);
    p = reshape(p, xs);
  else
    p = (x >= a & x <= b).*((ncdf(eta) - ncdf(alpha))/(ncdf(beta) - ncdf(alpha))).^n;
  end
else
  if length(m) > 1 || length(s) > 1
    tmpa = repmat(vpa_ncdf(alpha), size(x,1), 1);
    tmpb = repmat(vpa_ncdf(beta), size(x,1), 1);
    p = mean(((vpa_ncdf(eta) - tmpa)./(tmpb - tmpa)).^n, 2);
    p = reshape(p, xs);
  else
    p = (x >= a & x <= b).*((vpa_ncdf(eta) - vpa_ncdf(alpha))/(vpa_ncdf(beta) - vpa_ncdf(alpha))).^n;
  end
  p = min(max(p, 10^-digits), vpa(1)-10^-digits);
end

end

