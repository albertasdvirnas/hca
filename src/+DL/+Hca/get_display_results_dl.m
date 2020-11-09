function [] = get_display_results_dl(barcodeGen, consensusStruct,comparisonStruct,theoryStruct, sets)
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
markers = ['o';'s';'x';'+';'d';'v';'^';'<';'>'];

% compute cummulative sum of lengths of barcodes
lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));

% how many barcodes were compared to
numBar = length(comparisonStruct);

% If consensus was generated, then actual number of barcodes is one less
if sets.genConsensus == 1
  numBar = numBar-length(consensusStruct);
end

maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',false)');
maxcoefDense = cell2mat(cellfun(@(x) x.maxcoefDense,comparisonStruct,'UniformOutput',false)');
maxcoefSparse = cell2mat(cellfun(@(x) x.maxcoefSparse,comparisonStruct,'UniformOutput',false)');
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
plot_best_bar_dl(fig1,barcodeGen(1,:),consensusStruct,comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, ii, 'dual', 1);

subplot(2,2,4); hold on
%todo: improve this plot with more information
theoryStruct2 = theoryStruct;
for i=1:length(theoryStruct)
  theoryStruct2{i}.filename = theoryStruct{i}.filename2;
end
plot_best_bar_dl(fig1,barcodeGen(2,:),consensusStruct,comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, ii, 'dual', 2);


ax = gca;
set(gca, 'LooseInset', get(gca,'TightInset'))
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig1.PaperPositionMode = 'auto';
fig_pos = fig1.PaperPosition;
fig1.PaperSize = [fig_pos(3) fig_pos(4)];

%%

% for i=1:numBar
%   [maxAll, b] = max([maxcoef(i,1) maxcoefDense(i,1) maxcoefSparse(i,1)]);
%   if maxAll < 3; continue; end
%   disp([i maxAll b])
%   figure;
%   subplot(2,1,1); hold on
%   switch b
%     case 1
%       plot_best_bar_dl(fig1,barcodeGen(1,:),consensusStruct,comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'dual', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),consensusStruct,comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'dual', 2);
%     case 2
%       plot_best_bar_dl(fig1,barcodeGen(1,:),consensusStruct,comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'dense', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),consensusStruct,comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'dense', 2);
%     case 3
%       plot_best_bar_dl(fig1,barcodeGen(1,:),consensusStruct,comparisonStruct, theoryStruct, sets.userDefinedSeqCushion, i, 'sparse', 1);
%       subplot(2,1,2);
%       plot_best_bar_dl(fig1,barcodeGen(2,:),consensusStruct,comparisonStruct, theoryStruct2, sets.userDefinedSeqCushion, i, 'sparse', 2);
%   end
%   hold off
% end


%% additional things/options, added in HCA 4.1

% take a simple example for concentric plot
%     idx = 87;


% todo: make more user friendly..
%     try
%         mkdir(sets.output.matDirpath,sets.timestamp);
%     end
%
%     try
%         saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'result_plot.eps'),'epsc');
%     end
%
%     try
%         if sets.plotallmatches == 1
%              mkdir(fullfile(sets.output.matDirpath,sets.timestamp),'Plots');
%             for i=1:size(maxcoef,1)
%                 max2 = nan(size(maxcoef));      max2(i,1) = maxcoef(i,1);
%                 fig1 = figure('Visible', 'off');
% %                 fig1 = figure;
%                 if isequal(sets.comparisonMethod,'mp') || isequal(sets.comparisonMethod,'mpAll') || isequal(sets.comparisonMethod,'hmm')
%                     ax1 = subplot(1,1,1);
%                     if max2~=0
%                         plot_best_bar_mp(ax1,barcodeGen,[],comparisonStruct, theoryStruct, max2,0,sets);
%                     end
%                 else
%                     plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, max2);
%                 end
%                 % mp_based_on_output_pcc_test
%                 saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([num2str(i) '_plot.jpg'])));
%
% %                 saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([sets.timestamp '_' num2str(i) '_plot.eps'])),'epsc');
%
%             end
%         end
%     end
%    assignin('base','hcaSessionStruct',hcaSessionStruct)

%  cache('hcaSessionStruct') = hcaSessionStruct ;
%     end
end
