function [fig1,maxcoef] = plot_max_coef_dl( fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers )
% plot_max_coef - pltos three maximum coefficients

numToPlot = 2;
r = plot(maxcoefSparse(:, 1:numToPlot),1:size(maxcoefSparse,1),'ok', 'MarkerSize', 6);
hold on
q = plot(maxcoefDense(:, 1:numToPlot),1:size(maxcoefDense,1),'or', 'MarkerSize', 6);
p = plot(maxcoef(:, 1:numToPlot),1:size(maxcoef,1),'ob', 'MarkerSize', 8);
for i=1:size(p, 1)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i+3);
  r(i).Marker = markers(i+6);
end

ylabel('Barcode nr.','Interpreter','latex')
xlabel('Maximum match score','Interpreter','latex')
try
  xlim([ ...
    min([ ...
    maxcoef(:, 1:numToPlot); ...
    maxcoefDense(:, 1:numToPlot); ...
    maxcoefSparse(:, 1:numToPlot)], [], 'all') ...
    max([ ...
    maxcoef(:, 1:numToPlot); ...
    maxcoefDense(:, 1:numToPlot); ...
    maxcoefSparse(:, 1:numToPlot)], [], 'all')]);
catch
  xlim([0.5 1]);
end
ylim([0,size(maxcoef,1)+2])
legText = {'$\hat C^{dots}$','$C_2^{dots}$','$C_3^{dots}$', ...
  '$\hat C^{cb}$','$C_2^{cb}$','$C_3^{cb}$', ...
  '$\hat C^{dual}$','$C_2^{dual}$','$C_3^{dual}$'};
legend(legText(1:4-numToPlot:end),'Location','southwest','Interpreter','latex')
end

