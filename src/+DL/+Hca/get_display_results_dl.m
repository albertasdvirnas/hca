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
lengthBorders = cumsum(cellfun(@(x) x.length, theoryStruct));

% how many barcodes were compared to
numBar = length(comparisonStruct);

maxcoef = cell2mat(cellfun(@(x) x.dual.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefDense = cell2mat(cellfun(@(x) x.dense.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefSparse = cell2mat(cellfun(@(x) x.sparse.maxcoef,comparisonStruct,'UniformOutput',false)');
[~, ii] = max(maxcoef(:,1));

% plot best positions
[molTableFileName, dirpath] = uigetfile('*.*', 'Select file containing molecule data-table');
if not(isequal(dirpath, 0))
  molDataTable = readtable(fullfile(dirpath, molTableFileName));
  molIds = cellfun(@(x) str2double(x.name(1:end-4)), barcodeGen(1,:));
  [~, tableIds] = ismember(molIds, molDataTable.MoleculeID);
  bionanoConfidence = molDataTable.Confidence(tableIds);
  bionanoZ = -norminv(10.^-bionanoConfidence);
  bionanoPositionsBp = min(molDataTable.RefStartPos(tableIds), ...
    molDataTable.RefStartPos(tableIds));
  bionanoPositions = bionanoPositionsBp/sets.pvalue.pixelWidth_nm*sets.pvalue.nmbp;
  bionanoIdx = 5; % TEMP, what do if theoryStruct is not sorted??!

  % plot max corr coefs
  subplot(2,2,1);hold on;
  import DL.Hca.plot_max_coef_dl;
  fig1 = plot_max_coef_dl(fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers, bionanoZ);

  subplot(2,2,2);hold on;
  import DL.Hca.plot_best_pos_dl;
  plot_best_pos_dl(fig1,comparisonStruct, numBar, sets, markers,lengthBorders, bionanoPositions, bionanoIdx);
else

  % plot max corr coefs
  subplot(2,2,1);hold on;
  import DL.Hca.plot_max_coef_dl;
  fig1 = plot_max_coef_dl(fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers);

  subplot(2,2,2);hold on;
  import DL.Hca.plot_best_pos_dl;
  plot_best_pos_dl(fig1, comparisonStruct, numBar, sets, markers, lengthBorders);
end

subplot(2,2,3); hold on
%todo: improve this plot with more information
import DL.Hca.plot_best_bar_dl;
plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 1);

subplot(2,2,4); hold on
%todo: improve this plot with more information
plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 2);


%%
% b = 3;
% btype = {'dual', 'dense', 'sparse'};
% for i=1:numBar
% %   [maxAll, b] = max([maxcoef(i,1) maxcoefDense(i,1) maxcoefSparse(i,1)]);
% %   if maxAll < 3; continue; end
% %   disp([i maxAll b])
%   figure;
%   subplot(2,1,1); hold on
%   plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, btype{b}, 1);
%   subplot(2,1,2);
%   plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, btype{b}, 2);
%   hold off
% end
end
