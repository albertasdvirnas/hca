function [m_fit, ms_fit, s_fit, n_fit] = trunc_normal_ev_fit(y, a, b, lowerBound, upperBound, startParams, isConstant, numQuenched)
if nargin < 2
  a = 0;
end
if nargin < 3
  b = 1;
end
if nargin < 4
  lowerBound = [a 1e-16 1e-16 1];
end
if nargin < 5
  upperBound = [b 1 1 1e16];
end
if nargin < 6
  startParams = min(max([(b+a)/2 0.1 0.1 1], lowerBound), upperBound);
end
if nargin < 7
  isConstant = false(1,4);
end
if nargin < 8
  numQuenched = 100;
end

pM = linspace(0.001, 0.999, numQuenched);

import DL.Hca.trunc_normal_ev_pdf
import DL.Hca.trunc_normal_ev_cdf
switch bin2dec(num2str(isConstant))
  case 0
    pdfit = @(data, m, ms, s, n) trunc_normal_ev_pdf(data, norminv(pM, m, ms), s, n, a, b);
    cdfit = @(data, m, ms, s, n) trunc_normal_ev_cdf(data, norminv(pM, m, ms), s, n, a, b);
  case 1
    pdfit = @(data, m, ms, s) trunc_normal_ev_pdf(data, norminv(pM, m, ms), s, startParams(4), a, b);
    cdfit = @(data, m, ms, s) trunc_normal_ev_cdf(data, norminv(pM, m, ms), s, startParams(4), a, b);
  case 2
    pdfit = @(data, m, ms, n) trunc_normal_ev_pdf(data, norminv(pM, m, ms), startParams(3), n, a, b);
    cdfit = @(data, m, ms, n) trunc_normal_ev_cdf(data, norminv(pM, m, ms), startParams(3), n, a, b);
  case 3
    pdfit = @(data, m, ms) trunc_normal_ev_pdf(data, norminv(pM, m, ms), startParams(3), startParams(4), a, b);
    cdfit = @(data, m, ms) trunc_normal_ev_cdf(data, norminv(pM, m, ms), startParams(3), startParams(4), a, b);
  case 4
    pdfit = @(data, m, s, n) trunc_normal_ev_pdf(data, norminv(pM, m, startParams(2)), s, n, a, b);
    cdfit = @(data, m, s, n) trunc_normal_ev_cdf(data, norminv(pM, m, startParams(2)), s, n, a, b);
  case 5
    pdfit = @(data, m, s) trunc_normal_ev_pdf(data, norminv(pM, m, startParams(2)), s, startParams(4), a, b);
    cdfit = @(data, m, s) trunc_normal_ev_cdf(data, norminv(pM, m, startParams(2)), s, startParams(4), a, b);
  case 6
    pdfit = @(data, m, n) trunc_normal_ev_pdf(data, norminv(pM, m, startParams(2)), startParams(3), n, a, b);
    cdfit = @(data, m, n) trunc_normal_ev_cdf(data, norminv(pM, m, startParams(2)), startParams(3), n, a, b);
  case 7
    pdfit = @(data, m) trunc_normal_ev_pdf(data, norminv(pM, m, startParams(2)), startParams(3), startParams(4), a, b);
    cdfit = @(data, m) trunc_normal_ev_cdf(data, norminv(pM, m, startParams(2)), startParams(3), startParams(4), a, b);
  case 8
    pdfit = @(data, ms, s, n) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), ms), s, n, a, b);
    cdfit = @(data, ms, s, n) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), ms), s, n, a, b);
  case 9
    pdfit = @(data, ms, s) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), ms), s, startParams(4), a, b);
    cdfit = @(data, ms, s) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), ms), s, startParams(4), a, b);
  case 10
    pdfit = @(data, ms, n) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), ms), startParams(3), n, a, b);
    cdfit = @(data, ms, n) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), ms), startParams(3), n, a, b);
  case 11
    pdfit = @(data, ms) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), ms), startParams(3), startParams(4), a, b);
    cdfit = @(data, ms) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), ms), startParams(3), startParams(4), a, b);
  case 12
    pdfit = @(data, s, n) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), startParams(2)), s, n, a, b);
    cdfit = @(data, s, n) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), startParams(2)), s, n, a, b);
  case 13
    pdfit = @(data, s) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), startParams(2)), s, startParams(4), a, b);
    cdfit = @(data, s) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), startParams(2)), s, startParams(4), a, b);
  case 14
    pdfit = @(data, n) trunc_normal_ev_pdf(data, norminv(pM, startParams(1), startParams(2)), startParams(3), n, a, b);
    cdfit = @(data, n) trunc_normal_ev_cdf(data, norminv(pM, startParams(1), startParams(2)), startParams(3), n, a, b);
  case 15
    error("At least one parameter must not be constant.")
end

warning('off', 'stats:mlecov:NonPosDefHessian')
warning('off', 'stats:mle:IterLimit')
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'TolX',1e-6);%'Display','iter');
[params, ~] = mle(...
  y(~isnan(y)), 'pdf', pdfit, ...
  'start', startParams(not(isConstant)), ...
  'cdf', cdfit, ...
  'Lowerbound', lowerBound(not(isConstant)), ...
  'Upperbound', upperBound(not(isConstant)), ...
  'Alpha', 0.01, ...
  'Options', opt);
warning('on', 'stats:mlecov:NonPosDefHessian')
warning('on', 'stats:mle:IterLimit')

m_fit = startParams(1);
ms_fit = startParams(2);
s_fit = startParams(3);
n_fit = startParams(4);
switch bin2dec(num2str(isConstant))
  case 0
    m_fit = params(1);
    ms_fit = params(2);
    s_fit = params(3);
    n_fit = params(4);
  case 1
    m_fit = params(1);
    ms_fit = params(2);
    n_fit = params(3);
  case 2
    m_fit = params(1);
    ms_fit = params(2);
    s_fit = params(3);
  case 3
    m_fit = params(1);
    ms_fit = params(2);
  case 4
    m_fit = params(1);
    s_fit = params(2);
    n_fit = params(3);
  case 5
    m_fit = params(1);
    s_fit = params(2);
  case 6
    m_fit = params(1);
    n_fit = params(2);
  case 7
    m_fit = params(1);
  case 8
    ms_fit = params(1);
    s_fit = params(2);
    n_fit = params(3);
  case 9
    ms_fit = params(1);
    s_fit = params(2);
  case 10
    ms_fit = params(1);
    n_fit = params(2);
  case 11
    ms_fit = params(1);
  case 12
    s_fit = params(1);
    n_fit = params(2);
  case 13
    s_fit = params(1);
  case 14
    n_fit = params(1);
end
