function [fig1] = plot_best_pos_dl( fig1,comparisonStruct, numBar, sets, markers,lengthBorders)
% plot_best_pos - pltos three maximum coefficients

cumLengths = [0 lengthBorders]';

posShift = cumLengths(cellfun(@(x) x.idx,comparisonStruct)');

numToPlot = 2;
pos = cell2mat(cellfun(@(x) x.dual.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');
posDense = cell2mat(cellfun(@(x) x.dense.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');
posSparse = cell2mat(cellfun(@(x) x.sparse.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');

r = plot(posSparse+posShift,1:size(posSparse,1),'ok', 'MarkerSize', 6);
hold on
q = plot(posDense+posShift,1:size(posDense,1),'or', 'MarkerSize', 6);
p = plot(pos+posShift,1:size(pos,1),'ob', 'MarkerSize', 8);
for i=1:size(p, 1)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i+3);
  r(i).Marker = markers(i+6);
end
hold on

plot(lengthBorders,zeros(1,length(lengthBorders)),'x')

xlabel('Best position (px)','Interpreter','latex')
ylabel('Barcode nr.','Interpreter','latex')
legText = {'$\hat C^{dots}$','$C_2^{dots}$','$C_3^{dots}$', ...
  '$\hat C^{cb}$','$C_2^{cb}$','$C_3^{cb}$', ...
  '$\hat C^{dual}$','$C_2^{dual}$','$C_3^{dual}$'};
legend(legText(1:4-numToPlot:end),'Location','southwest','Interpreter','latex')

ylim([0,size(pos,1)+2])

end



