function [fig1,maxcoef] = plot_max_coef_dl( fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers, bionanoCoef )
% plot_max_coef - pltos three maximum coefficients

numToPlot = 2;
r = plot(maxcoefSparse(:, 1:numToPlot),1:size(maxcoefSparse,1),'ok', 'MarkerSize', 6);
hold on
q = plot(maxcoefDense(:, 1:numToPlot),1:size(maxcoefDense,1),'or', 'MarkerSize', 6);
p = plot(maxcoef(:, 1:numToPlot),1:size(maxcoef,1),'ob', 'MarkerSize', 8);
if nargin > 7
  plot(bionanoCoef, 1:size(maxcoef,1),'.','Color', [0 .8 .5], 'MarkerSize', 10)
end
for i=1:size(p, 1)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i+3);
  r(i).Marker = markers(i+6);
end

ylabel('Barcode nr.','Interpreter','latex')
xlabel('Maximum match score','Interpreter','latex')
try
  xbounds = [ ...
    min([ ...
    maxcoef(:, 1:numToPlot); ...
    maxcoefDense(:, 1:numToPlot); ...
    maxcoefSparse(:, 1:numToPlot)], [], 'all') ...
    max([ ...
    maxcoef(:, 1:numToPlot); ...
    maxcoefDense(:, 1:numToPlot); ...
    maxcoefSparse(:, 1:numToPlot)], [], 'all')];
  if nargin > 7
    xbounds(1) = min(xbounds(1), min(bionanoCoef));
    xbounds(2) = max(xbounds(2), max(bionanoCoef));
  end
  xlim(xbounds);
catch
  xlim([0.5 1]);
end
ylim([0,size(maxcoef,1)+2])
legText = [ ...
  arrayfun(@(i) strcat("$\hat Z^{dots}_", num2str(i), "$"), 1:numToPlot, 'un', 0) ...
  arrayfun(@(i) strcat("$\hat Z^{cb}_", num2str(i), "$"), 1:numToPlot, 'un', 0) ...
  arrayfun(@(i) strcat("$\hat Z^{dual}_", num2str(i), "$"), 1:numToPlot, 'un', 0)
  ];
if nargin > 7
  legText{end+1} = 'Bionano positions';
end
legend(legText,'Location','southwest','Interpreter','latex')
end

