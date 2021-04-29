function [] = get_display_results_dl(barcodeGen, comparisonStruct,theoryStruct, sets, externalAlignmentStruct)
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
maxcoefDense = cell2mat(cellfun(@(x) x.ch1.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefSparse = cell2mat(cellfun(@(x) x.ch2.maxcoef,comparisonStruct,'UniformOutput',false)');
[~, ii] = max(maxcoef(:,1));

% plot best positions

if nargin > 5 && not(isempty(fieldnames(externalAlignmentStruct)))
  % plot max corr coefs
  subplot(2,2,1);hold on;
  import DL.Hca.plot_max_coef_dl;
  fig1 = plot_max_coef_dl(fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers, 1, externalAlignmentStruct.Z);
  
  subplot(2,2,2);hold on;
  import DL.Hca.plot_best_pos_dl;
  plot_best_pos_dl(fig1,comparisonStruct, numBar, sets, markers,lengthBorders, 1, externalAlignmentStruct.pos, externalAlignmentStruct.idx);
else
  
  % plot max corr coefs
  subplot(2,2,1);hold on;
  import DL.Hca.plot_max_coef_dl;
  fig1 = plot_max_coef_dl(fig1, maxcoef, maxcoefDense, maxcoefSparse, numBar, sets, markers, 1);
  
  subplot(2,2,2);hold on;
  import DL.Hca.plot_best_pos_dl;
  plot_best_pos_dl(fig1, comparisonStruct, numBar, sets, markers, lengthBorders, 1);
end

subplot(2,2,3); hold on
%todo: improve this plot with more information
import DL.Hca.plot_best_bar_dl;
plot_best_bar_dl(fig1,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 1);

subplot(2,2,4); hold on
%todo: improve this plot with more information
plot_best_bar_dl(fig1,barcodeGen(2,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 2);
hold off

%%
% %
% stretchF = cellfun(@(x) x.external.bestBarStretch, comparisonStruct);
% for i=1:numBar
%   import DL.Hca.plot_bar_at_pos;
%   f = figure;
%   subplot(2,1,1); hold on
%   if isfield(comparisonStruct{i}.external, 'maxcoefParts')
%     plot_bar_at_pos(barcodeGen(1,:), i, theoryStruct, externalAlignmentStruct.pos(i), externalAlignmentStruct.or(i), stretchF(i), externalAlignmentIdx, comparisonStruct{i}.external.maxcoefPartsCC(1), sets.userDefinedSeqCushion, 1)
%     subplot(2,1,2);
%     plot_bar_at_pos(barcodeGen(2,:), i, theoryStruct, externalAlignmentStruct.pos(i), externalAlignmentStruct.or(i), stretchF(i), externalAlignmentIdx, comparisonStruct{i}.external.maxcoefPartsCC(2), sets.userDefinedSeqCushion, 2)
%   else
%     plot_bar_at_pos(barcodeGen(1,:), i, theoryStruct, externalAlignmentStruct.pos(i), externalAlignmentStruct.or(i), stretchF(i), externalAlignmentIdx, nan, sets.userDefinedSeqCushion, 1)
%     subplot(2,1,2);
%     plot_bar_at_pos(barcodeGen(2,:), i, theoryStruct, externalAlignmentStruct.pos(i), externalAlignmentStruct.or(i), stretchF(i), externalAlignmentIdx, nan, sets.userDefinedSeqCushion, 2)
%   end
%   hold off
% end



% b = 1;
% btype = {'dual', 'ch1', 'ch2'};
% for i=1:numBar
% %   [maxAll, b] = max([maxcoef(i,1) maxcoefDense(i,1) maxcoefSparse(i,1)]);
% %   if maxcoefDense(i,1) < 3; continue; end
%   disp(i)
% %   disp([i maxAll b])
%   f = figure;
%   subplot(2,1,1); hold on
%   plot_best_bar_dl(f,barcodeGen(1,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, btype{b}, 1);
%   subplot(2,1,2);
%   plot_best_bar_dl(f,barcodeGen(2,:),comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, btype{b}, 2);
%   hold off
%   saveas(f, fullfile(pwd, 'out', barcodeGen{1, i}.name(1:end-4)), 'png');
%   close(f)
% end
end
