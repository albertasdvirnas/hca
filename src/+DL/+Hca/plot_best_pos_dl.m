function [fig1] = plot_best_pos_dl( fig1,comparisonStruct, numBar, sets, markers,lengthBorders, bionanoPositions, bionanoIdx)
% plot_best_pos - pltos three maximum coefficients

cumLengths = [0 lengthBorders]';

posShift = cumLengths(cellfun(@(x) x.dual.idx,comparisonStruct)');
posShiftDense = cumLengths(cellfun(@(x) x.dense.idx,comparisonStruct)');
posShiftSparse = cumLengths(cellfun(@(x) x.sparse.idx,comparisonStruct)');

numToPlot = 2;
pos = cell2mat(cellfun(@(x) x.dual.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');
posDense = cell2mat(cellfun(@(x) x.dense.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');
posSparse = cell2mat(cellfun(@(x) x.sparse.pos(1:numToPlot),comparisonStruct,'UniformOutput',0)');

r = plot(posSparse+posShiftSparse,1:size(posSparse,1),'ok', 'MarkerSize', 6);
hold on
q = plot(posDense+posShiftDense,1:size(posDense,1),'or', 'MarkerSize', 6);
p = plot(pos+posShift,1:size(pos,1),'ob', 'MarkerSize', 8);
if nargin > 6
  posShiftBionano = cumLengths(bionanoIdx);
  plot(bionanoPositions+posShiftBionano, 1:size(pos,1),'.','Color', [0 .8 .5], 'MarkerSize', 10)
end
for i=1:size(p, 1)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i+3);
  r(i).Marker = markers(i+6);
end
hold on

plot(lengthBorders,zeros(1,length(lengthBorders)),'x')

xlabel('Best position (px)','Interpreter','latex')
ylabel('Barcode nr.','Interpreter','latex')
legText = [ ...
  arrayfun(@(i) strcat("$\hat Z^{dots}_", num2str(i), "$"), 1:numToPlot, 'un', 0) ...
  arrayfun(@(i) strcat("$\hat Z^{cb}_", num2str(i), "$"), 1:numToPlot, 'un', 0) ...
  arrayfun(@(i) strcat("$\hat Z^{dual}_", num2str(i), "$"), 1:numToPlot, 'un', 0)
  ];
if nargin > 6
  legText{end+1} = 'Bionano positions';
end
legend(legText,'Location','northeast','Interpreter','latex')
ylim([0,size(pos,1)+2])

end



