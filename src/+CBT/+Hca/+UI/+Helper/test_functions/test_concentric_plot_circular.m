load('concentric_comparison_test.mat');

maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',false)');

    
theoryStruct{1}.isLinearTF = 0;

sets.theory.isLinearTF = theoryStruct{1}.isLinearTF;
fig=figure;
hAxis = subplot(1,1,1);
import CBT.Hca.UI.Helper.plot_best_concentric_image;
plot_best_concentric_image(hAxis,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets);
