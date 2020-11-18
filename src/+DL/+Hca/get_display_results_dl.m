function [] = get_display_results_dl(barcodeGen, comparisonStruct,theoryStruct, sets)
% get_display_results
% Display results from the comparison of experiments vs theory
%     Args:
%         barcodeGen: Input
%         consensusStruct:
%         comparisonStruct:
%         theoryStruct:
%         sets:

% if it was chosen to display results
if sets.displayResults ==0
  fig1 = figure('Visible', 'off');
else
  fig1 = figure;
end

% choose markers for everything
markers = ['o';'s';'x';'v';'+';'d';'^';'<';'>'];

% compute cummulative sum of lengths of barcodes
lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));

% how many barcodes were compared to
numBar = length(comparisonStruct);

maxcoef = cell2mat(cellfun(@(x) x.dual.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefDense = cell2mat(cellfun(@(x) x.dense.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefSparse = cell2mat(cellfun(@(x) x.sparse.maxcoef,comparisonStruct,'UniformOutput',false)');
[~, ii] = max(maxcoef(:,1));

% plot max corr coefs
subplot(2,2,1);hold on;
import DL.Hca.plot_max_coef_dl;
fig1 = plot_max_coef_dl(fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers);

% plot best positions
subplot(2,2,2);hold on;
import DL.Hca.plot_best_pos_dl;
plot_best_pos_dl(fig1,comparisonStruct, numBar, sets, markers,lengthBorders);

subplot(2,2,3); hold on
%todo: improve this plot with more information
import DL.Hca.plot_best_bar_dl;
plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 1);

subplot(2,2,4); hold on
%todo: improve this plot with more information
theoryStruct2 = theoryStruct;
for i=1:length(theoryStruct)
  theoryStruct2{i}.filename = theoryStruct{i}.filename2;
end
plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, ii, 'dual', 2);


%%

% for i=1:numBar
%   [maxAll, b] = max([maxcoef(i,1) maxcoefDense(i,1) maxcoefSparse(i,1)]);
% %   if maxAll < 3; continue; end
%   disp([i maxAll b])
%   figure;
%   subplot(2,1,1); hold on
%   switch b
%     case 1
%       plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'dual', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'dual', 2);
%     case 2
%       plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'dense', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'dense', 2);
%     case 3
%       plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'sparse', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'sparse', 2);
%   end
%   hold off
% end
end
