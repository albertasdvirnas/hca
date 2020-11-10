function [p] = beta_ev_cdf(x, a, b, n, extraPrecision)

if length(a) > 1
  xs = size(x);
  x = x(:);
  a = a(:)';
end

if nargin < 5 || ~extraPrecision
  if length(a) > 1
    p = mean(cell2mat(arrayfun(@(a1) (0.5*(1 + betainc(x.^2, 0.5*b, a1/2-1, 'lower'))).^n, a, 'un', 0)), 2);
    p = reshape(p, xs);
  else
    p = (0.5*(1 + betainc(x.^2, 0.5*b, a/2-1, 'lower'))).^n;
  end
  p = min(max(p, 10^-digits), 1-10^-digits);
else
  if length(a) > 1
    p = mean(cell2mat(arrayfun(@(a1) (0.5*(vpa(1) + vpa(betainc(x.^2, 0.5*b, a1/2-1, 'lower')))).^n, a, 'un', 0)), 2);
    p = reshape(p, xs);
  else
    p = (0.5*(vpa(1) + vpa(betainc(x.^2, 0.5*b, a/2-1, 'lower')))).^n;
  end
  p = min(max(p, 10^-digits), vpa(1)-10^-digits);
end

