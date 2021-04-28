function [p] = beta_ev_cdf(x, a, b, n, extraPrecision)

  undefx = x.^2 <= 0 | x.^2 >= 1 | isnan(x);
  x(undefx) = -1;

  if nargin < 5 ||~extraPrecision
    p = (0.5 * (1 + sign(x) .* betainc(x.^2, 0.5 * b, a / 2 - 1, 'lower'))).^n;
  else
    p = (0.5 * (vpa(2) - sign(x) .* betainc(x.^2, 0.5 * b, a / 2 - 1, 'upper'))).^n;
    p = min(p, vpa(1) - 10^ - digits);
  end

  p(undefx) = nan;

  if nargin > 4 && extraPrecision
    p(x == 0) = 10^ - digits;
    p(x == 1) = vpa(1) - 10^ - digits;
  else
    p(x == 0) = 0;
    p(x == 1) = 1;
  end
