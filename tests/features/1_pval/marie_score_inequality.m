% here we show how to apply Marie's method to HCA data

% first generate some data / these settings should be maybe set for this
% individually
% sets =ini2struct('sp_sim_sets_1.txt');
% setstxt = 'sp_settings.txt';

% % chech if the sequences have to be reproducible
rng(1,'twister');
    
% one kymo and one theory
numRows = 1;
numKymos = 2;
numTheories = 1;
totalLength = 10000; % total length human DNA chromosomes.
sets.lengthN = repmat(totalLength,1,numTheories);
sets.outFold = 'test/';
sets.strF = 2;
sets.pccScore = 0.9;
sets.lambda = 500;
sets.kernelsigma = 2.3;
[randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_multi_sim_theory_multi(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows,numKymos,numTheories);

% we interpolate to so that we callculate PCC every 100 basepairs (1/5th of
% a pixel, a distance at which theoretical barcode changes (approximately)
% ).
% 
% 
% vec = [1:length(randKymos{1})]';                    % Create Data
% xi = linspace(1,length(vec),length(vec)*5);           % Interpolation Vector
% data2 = interp1( vec , randKymos{1}, xi ,'spline') ;

% now take 125 kb for a barcode = 1250 points ( so rather long!)

% now run comparison

% now run comparison
numFrames =1;
nmBp = 0.3;
comparisonMethod = 'mp';
[hcaStruct] = run_parallel('hca_p1.txt', theoriesTxt, 'test/', numFrames, nmBp,comparisonMethod);

comparisonMethod = 'mass_pcc';
[hcaStructMASS] = run_parallel('hca_p1.txt', theoriesTxt, 'test/', numFrames, nmBp,comparisonMethod);

% extract single barcode based on run_parallel results (so we can quickly
% compute to find max corr coefficient) 

import generate.kymo_to_barcode;
barcodeGenNew = kymo_to_barcode(hcaStruct.barcodeGenC,cellfun(@(x) x.secondPos,hcaStruct.comparisonStruct),...
    cellfun(@(x) x.lengthMatch,hcaStruct.comparisonStruct),cellfun(@(x) x.bestBarStretch,hcaStruct.comparisonStruct));
% now convert this to m

% now run comparison for single barcode using mass_pcc (and outputing all
% the score, because we use these for Marie's distribution fit
comparisonMethod = 'mass_pcc';
[hcaStructN] = run_parallel_from_barcodes('hca_p1.txt', theoriesTxt,barcodeGenNew{1}, numFrames, nmBp,comparisonMethod);

% these have to be the same for all tests for method to work
hcaStructN.comparisonStruct{1}.maxcoef
hcaStruct.comparisonStruct{1}.maxcoef

% hcaStructN.comparisonStruct{1}.pos(1)-hcaStruct.comparisonStruct{1}.secondPos
% hcaStruct.comparisonStruct{1}.pos

data=hcaStructN.comparisonStruct{1}.dist(:);
data2 = sort(data(:));
pd = fitdist(data2(1:end-3),'normal');

x_values = -0.7:0.01:0.7;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
hold on
dat = histogram(data(:),'Normalization','pdf','DisplayStyle','Stairs')
set(gca, 'YScale', 'log')

xvals = (dat.BinEdges(1:end-1)+dat.BinEdges(2:end))/2;
yvals = dat.Values;

g = fittype('a*exp(-((x-b)/c)^2)+d*x^2');
FO = fit(xvals', yvals', g);

figure,
x_values = -0.7:0.01:0.7;
y = FO(x_values);
plot(x_values,y,'LineWidth',2)
hold on
dat = histogram(data(:),'Normalization','pdf','DisplayStyle','Stairs')
set(gca, 'YScale', 'log')



comparisonMethod = 'mpAll';
% comparisonMethod = 'mp';
% comparisonMethod = 'spearman';
% comparisonMethod = 'dtw';
[rezMax,bestBarStretch,bestLength,rezMaxAll] = run_parallel_str('hca_p1.txt', theoriesTxt, sets.outFold, numFrames, nmBp,comparisonMethod);

