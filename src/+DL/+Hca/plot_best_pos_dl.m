function [fig1] = plot_best_pos_dl( fig1,comparisonStruct, numBar, sets, markers,lengthBorders)
% plot_best_pos - pltos three maximum coefficients

cumLengths = [0 lengthBorders]';

posShift = cumLengths(cellfun(@(x) x.idx,comparisonStruct)');

pos = cell2mat(cellfun(@(x) x.dual.pos,comparisonStruct,'UniformOutput',0)');
posDense = cell2mat(cellfun(@(x) x.dense.pos,comparisonStruct,'UniformOutput',0)');
posSparse = cell2mat(cellfun(@(x) x.sparse.pos,comparisonStruct,'UniformOutput',0)');

r = plot(posSparse+posShift,1:size(posSparse,1),'ok', 'MarkerSize', 6);
hold on
q = plot(posDense+posShift,1:size(posDense,1),'or', 'MarkerSize', 6);
p = plot(pos+posShift,1:size(pos,1),'ob', 'MarkerSize', 8);
for i=1:length(p)
  p(i).Marker = markers(i);
  q(i).Marker = markers(i+3);
  r(i).Marker = markers(i+6);
end
hold on

% if  sets.genConsensus == 1
%   plot(0:100:sum(lengthBorders), 0.5+repmat(numBar,length(0:100:sum(lengthBorders)),1))
% end

plot(lengthBorders,zeros(1,length(lengthBorders)),'x')

xlabel('Best position (px)','Interpreter','latex')
ylabel('Barcode nr.','Interpreter','latex')
if size(pos,2) == 1
  legend({'$\hat C$','Theories seperator'},'Location','southwest','Interpreter','latex')
else
  legend({'$\hat C^{dots}$','$C_2^{dots}$','$C_3^{dots}$', ...
    '$\hat C^{cb}$','$C_2^{cb}$','$C_3^{cb}$', ...
    '$\hat C^{dual}$','$C_2^{dual}$','$C_3^{dual}$'},'Location','southwest','Interpreter','latex')
end

ylim([0,size(pos,1)+2])

end



