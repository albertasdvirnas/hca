function p = trunc_normal_ev_cdf(x, m, s, n, a, b, extraPrecision)
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
alph = (a - m)./s;
bet = (b - m)./s;

ncdf = @(x) .5*(1+erf(x/sqrt(2)));
vpa_ncdf = @(x) .5*(vpa(2)-erfc(x/sqrt(2)));
alph2 = ncdf(alph);
bet2 = ncdf(bet);

if nargin < 7 || ~extraPrecision
  if length(m) > 1 || length(s) > 1
    p = mean(((ncdf(eta) - alph2)./(bet2 - alph2)).^n, 2);
    p = reshape(p, xs);
  else
    p = ((ncdf(eta) - alph2)/(bet2 - alph2)).^n;
  end
else
  if length(m) > 1 || length(s) > 1
    % This operation is extremely slow using vpa :(
    bet3 = repmat(bet2 - alph2, length(x), 1);
    alph3 = alph2./bet3;
    p = mean((vpa_ncdf(eta)./bet3 - alph3).^n, 2);
    p = reshape(p, xs);
  else
    p = ((vpa_ncdf(eta) - alph2)/(bet2 - alph2)).^n;
  end
  p = min(p, vpa(1)-10^-digits);
end