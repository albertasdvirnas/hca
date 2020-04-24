load('concentric_comparison_test.mat');

maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',false)');

    
theoryStruct{1}.isLinearTF = 0;


sets.output.userDefinedSeqCushion = 20;

fig1 = figure;
subplot(1,1,1), hold on
import CBT.Hca.UI.Helper.plot_best_bar;
plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets.output.userDefinedSeqCushion);

mkdir(fullfile(sets.output.matDirpath,sets.timestamp),'Plots');
for i=1:size(maxcoef,1)
    max2 = nan(size(maxcoef));      max2(i,1) = maxcoef(i,1);
    fig1 = figure('Visible', 'off');
    plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, max2);
    saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([num2str(i) '_plot.jpg'])));
end
            
  
%     
% fig=figure;
% hAxis = subplot(1,1,1);
% import CBT.Hca.UI.Helper.plot_best_concentric_image;
% plot_best_concentric_image(hAxis,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets);
