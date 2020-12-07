function [p] = beta_ev_cdf(x, a, b, n, extraPrecision)

undefx = x.^2 <= 0 | x.^2 >= 1 | isnan(x);
x(undefx) = 0.5;

if length(a) > 1
  xs = size(x);
  x = x(:);
  a = a(:)';
end

if nargin < 5 || ~extraPrecision
  if length(a) > 1
    p = mean(cell2mat(arrayfun(@(a1) (0.5*(1 + sign(x).*betainc(x.^2, 0.5*b, a1/2-1, 'lower'))).^n, a, 'un', 0)), 2);
    p = reshape(p, xs);
  else
    p = (0.5*(1 + sign(x).*betainc(x.^2, 0.5*b, a/2-1, 'lower'))).^n;
  end
else
  if length(a) > 1
    p = mean(cell2mat(arrayfun(@(a1) (0.5*(vpa(2) - sign(x).*betainc(x.^2, 0.5*b, a1/2-1, 'upper'))).^n, a, 'un', 0)), 2);
    p = reshape(p, xs);
  else
    p = (0.5*(vpa(2) - sign(x).*betainc(x.^2, 0.5*b, a/2-1, 'upper'))).^n;
  end
  p = min(p, vpa(1)-10^-digits);
end

p(undefx) = nan;