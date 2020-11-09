function [fig1,maxcoef] = plot_max_coef_dl( fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers )
% plot_max_coef - pltos three maximum coefficients

p = plot(maxcoef,1:size(maxcoef,1),'ob', 'MarkerSize', 8);
hold on
q = plot(maxcoefDense,1:size(maxcoefDense,1),'or', 'MarkerSize', 6);
r = plot(maxcoefSparse,1:size(maxcoefSparse,1),'ok', 'MarkerSize', 4);
for i=1:length(p)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i);
  r(i).Marker = markers(i);
end


% if  sets.genConsensus == 1
%   plot(0.1:0.1:1, 0.5+repmat(numBar,10,1));
% end

ylabel('Barcode nr.','Interpreter','latex')
xlabel('Maximum match score','Interpreter','latex')
try
  xlim([ ...
    min([ ...
    maxcoef(:); ...
    maxcoefDense(:); ...
    maxcoefSparse(:)]) ...
    max([ ...
    maxcoef(:); ...
    maxcoefDense(:); ...
    maxcoefSparse(:)])]);
catch
  xlim([0.5 1]);
end
ylim([0,size(maxcoef,1)+2])
if size(maxcoef,2) == 1
  legend({'$\hat C$','Theories separator'},'Location','southwest','Interpreter','latex')
else
  legend({'$\hat C^{dual}$','$C_2^{dual}$','$C_3^{dual}$', ...
    '$\hat C^{cb}$','$C_2^{cb}$','$C_3^{cb}$', ...
    '$\hat C^{dots}$','$C_2^{dots}$','$C_3^{dots}$'},'Location','southwest','Interpreter','latex')
end
end

