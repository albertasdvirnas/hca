function [p] = beta_ev_pdf(x, a, b, n)

if length(a) > 1
  x = x(:);
  a = a(:)';
end

if length(a) > 1
  p = mean(cell2mat(arrayfun(@(a1) n.*((.5*(1 + betainc(x.^2, .5*b, a1/2 - 1, 'lower'))).^(n-1)).*x.^(b-1).*((1-x.^2).^(a1/2-2))./beta(.5*b, a1/2-1), a, 'un', 0)), 2);
else
  p = n.*((.5*(1 + betainc(x.^2, .5*b, a/2 - 1, 'lower'))).^(n-1)).*x.^(b-1).*((1-x.^2).^(a/2-2))/beta(.5*b, a/2-1);
end
p = max(p, 10^-digits);

