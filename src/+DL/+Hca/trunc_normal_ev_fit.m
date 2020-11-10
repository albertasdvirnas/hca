function [m_fit, n_fit, s_fit] = trunc_normal_ev_fit(y, a, b, constantVal, lowerBound, upperBound, startParams)
if nargin < 2
  a = -1;
end
if nargin < 3
  b = 1;
end
if nargin < 4
  constantVal = nan(1,3);
end
if nargin < 5
  lowerBound = [max(a, -1e16) 1 1e-16];
end
if nargin < 6
  upperBound = [min(b, 1e16) 1e+16 1e+16];
end
if nargin < 7
  startParams = min(max([0.5 2 1e-6], lowerBound), upperBound);
end
isConstant = not(isnan(constantVal));

if isConstant(1)
  if isConstant(3)
    pdfit = @(data, n, s) trunc_normal_ev_pdf(data, constantVal(1), n, constantVal(3), a, b);
    cdfit = @(data, n, s) trunc_normal_ev_cdf(data, constantVal(1), n, constantVal(3), a, b);
  elseif isConstant(2)
    pdfit = @(data, s) trunc_normal_ev_pdf(data, constantVal(1), constantVal(2), s, a, b);
    cdfit = @(data, s) trunc_normal_ev_cdf(data, constantVal(1), constantVal(2), s, a, b);
  else
    pdfit = @(data, n, s) trunc_normal_ev_pdf(data, constantVal(1), n, s, a, b);
    cdfit = @(data, n, s) trunc_normal_ev_cdf(data, constantVal(1), n, s, a, b);
  end
elseif isConstant(2)
  if isConstant(3)
    pdfit = @(data, m) trunc_normal_ev_pdf(data, m, constantVal(2), constantVal(3), a, b);
    cdfit = @(data, m) trunc_normal_ev_cdf(data, m, constantVal(2), constantVal(3), a, b);
  else
    pdfit = @(data, m, s) trunc_normal_ev_pdf(data, m, constantVal(2), s, a, b);
    cdfit = @(data, m, s) trunc_normal_ev_cdf(data, m, constantVal(2), s, a, b);
  end
elseif isConstant(3)
  pdfit = @(data, m, n) trunc_normal_ev_pdf(data, m, n, constantVal(3), a, b);
  cdfit = @(data, m, n) trunc_normal_ev_cdf(data, m, n, constantVal(3), a, b);
else
  pdfit = @(data, m, n, s) trunc_normal_ev_pdf(data, m, n, s, a, b);
  cdfit = @(data, m, n, s) trunc_normal_ev_cdf(data, m, n, s, a, b);
end


% warning('off', 'stats:mlecov:NonPosDefHessian')
%     warning('off', 'stats:mle:IterLimit')
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'TolX',1e-6);%'Display','iter');
[params, ~] = mle(...
  y(~isnan(y)), 'pdf', pdfit, ...
  'start', startParams(not(isConstant)), ...
  'cdf', cdfit, ...
  'Lowerbound', lowerBound(not(isConstant)), ...
  'Upperbound', upperBound(not(isConstant)), ...
  'Alpha', 0.01, ...
  'Options', opt);
% warning('on', 'stats:mlecov:NonPosDefHessian')
%     warning('on', 'stats:mle:IterLimit')

if isConstant(1)
  m_fit = constantVal(1);
  if isConstant(3)
    n_fit = params(1);
    s_fit = constantVal(3);
  elseif isConstant(2)
    n_fit = constantVal(2);
    s_fit = params(1);
  else
    n_fit = params(1);
    s_fit = params(2);
  end
else
  m_fit = params(1);
  if isConstant(2)
    if isConstant(3)
      n_fit = constantVal(2);
      s_fit = constantVal(3);
    else
      n_fit = constantVal(2);
      s_fit = params(2);
    end
  elseif isConstant(3)
    n_fit = params(2);
    s_fit = constantVal(3);
  else
    n_fit = params(2);
    s_fit = params(3);
  end
end
end

