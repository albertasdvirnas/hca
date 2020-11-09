function p = qGauss_ev_cdf(x, r, N, extraPrecision)

if length(r) > 1
  xs = size(x);
  x = x(:);
  r = r(:)';
end

if nargin < 4 || ~extraPrecision
  if length(r) > 1
    p = mean((.5*(1+erf(sqrt(2./r-1).*erfinv(x)))).^N, 2);
    p = reshape(p, xs);
  else
    p = (.5*(1+erf(sqrt(2/r-1)*erfinv(x)))).^N;
  end
else
  if length(r) > 1
    p = mean((.5*(vpa(1)+erf(sqrt(vpa(2./r)-vpa(1)).*erfinv(x)))).^N, 2);
    p = reshape(p, xs);
  else
    p = (.5*(vpa(1)+erf(sqrt(vpa(2/r)-vpa(1))*erfinv(x)))).^N;
  end
  p = min(max(p, 10^-digits), vpa(1)-10^-digits);
end