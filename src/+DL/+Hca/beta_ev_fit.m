function [a_fit, b_fit, n_fit] = beta_ev_fit(y, lowerBound, upperBound, startParams, isConstant)

if nargin < 2
    lowerBound = [4 10^-digits 1];
end
if nargin < 3
    upperBound = [inf inf inf];
end
if nargin < 4
  startParams = min(max([4, 1, 1], lowerBound), upperBound);
end
if nargin < 5
    isConstant = false(1,3);
end

import DL.Hca.beta_ev_pdf
import DL.Hca.beta_ev_cdf
switch bin2dec(num2str(isConstant))
  case 0
    pdfit = @(data, a, b, n) beta_ev_pdf(data, a, b, n);
    cdfit = @(data, a, b, n) beta_ev_cdf(data, a, b, n);
  case 1
    pdfit = @(data, a, b) beta_ev_pdf(data, a, b, startParams(3));
    cdfit = @(data, a, b) beta_ev_cdf(data, a, b, startParams(3));
  case 2
    pdfit = @(data, a, n) beta_ev_pdf(data, a, startParams(2), n);
    cdfit = @(data, a, n) beta_ev_cdf(data, a, startParams(2), n);
  case 3
    pdfit = @(data, a) beta_ev_pdf(data, a, startParams(2), startParams(3));
    cdfit = @(data, a) beta_ev_cdf(data, a, startParams(2), startParams(3));
  case 4
    pdfit = @(data, b, n) beta_ev_pdf(data, startParams(1), b, n);
    cdfit = @(data, b, n) beta_ev_cdf(data, startParams(1), b, n);
  case 5
    pdfit = @(data, b) beta_ev_pdf(data, startParams(1), b, startParams(3));
    cdfit = @(data, b) beta_ev_cdf(data, startParams(1), b, startParams(3));
  case 6
    pdfit = @(data, n) beta_ev_pdf(data, startParams(1), startParams(2), n);
    cdfit = @(data, n) beta_ev_cdf(data, startParams(1), startParams(2), n);
  case 7
    error("At least one parameter must not be constant.")
end

warning('off', 'stats:mlecov:NonPosDefHessian')
warning('off', 'stats:mle:IterLimit')
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5, 'TolX', 1e-6);%,'Display','iter');
[params, ~] = mle(y(~isnan(y)), 'pdf', pdfit, ...
    'start', startParams(not(isConstant)), ...
    'cdf', cdfit, ...
    'Lowerbound', lowerBound(not(isConstant)), ...
    'Upperbound', upperBound(not(isConstant)), ...
    'Alpha', 0.01, ...
    'Options', opt);
warning('on', 'stats:mlecov:NonPosDefHessian')
warning('on', 'stats:mle:IterLimit')

a_fit = startParams(1);
b_fit = startParams(2);
n_fit = startParams(3);
switch bin2dec(num2str(isConstant))
  case 0
    a_fit = params(1);
    b_fit = params(2);
    n_fit = params(3);
  case 1
    a_fit = params(1);
    b_fit = params(2);
  case 2
    a_fit = params(1);
    n_fit = params(2);
  case 3
    a_fit = params(1);
  case 4
    b_fit = params(1);
    n_fit = params(2);
  case 5
    b_fit = params(1);
  case 6
    n_fit = params(1);
end
