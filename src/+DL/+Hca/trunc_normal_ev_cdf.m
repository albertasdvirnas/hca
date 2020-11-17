function p = trunc_normal_ev_cdf(x, m, s, n, a, b, extraPrecision)
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

ncdf = @(x) .5*(1+erf(x/sqrt(2)));
vpa_ncdf = @(x) .5*(vpa(1)+erf(vpa(x/sqrt(2))));

if nargin < 7 || ~extraPrecision
  if length(m) > 1 || length(s) > 1
    p = mean(((ncdf(eta) - ncdf(alpha))./(ncdf(beta) - ncdf(alpha))).^n, 2);
    p = reshape(p, xs);
  else
    p = ((ncdf(eta) - ncdf(alpha))/(ncdf(beta) - ncdf(alpha))).^n;
  end
else
  if length(m) > 1 || length(s) > 1
    % This operation is extremely slow using vpa :(
    tmpA = repmat(vpa_ncdf(alpha), xs(1), 1);
    tmpB = repmat(vpa_ncdf(beta), xs(1), 1);
    p = mean(((vpa_ncdf(eta) - tmpA)./(tmpB - tmpA)).^n, 2);
    p = reshape(p, xs);
  else
    p = ((vpa_ncdf(eta) - vpa_ncdf(alpha))/(vpa_ncdf(beta) - vpa_ncdf(alpha))).^n;
  end
  p = min(p, vpa(1)-10^-digits);
end

