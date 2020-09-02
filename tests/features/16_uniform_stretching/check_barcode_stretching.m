% this code is to check whether barcode is uniformally stretched. 

% first generate some data

sets =ini2struct('sp_sim_sets_1.txt');
setstxt = 'sp_settings.txt';

% % chech if the sequences have to be reproducible
rng(1,'twister');
    
% one kymo and one theory
numRows = 1;
numKymos = 2;
numTheories = 1;
sets.lengthN = repmat(10000,1,numTheories);
sets.outFold = 'test/';
sets.strF = 2;
sets.pccScore = 0.9;
sets.lambda = 1000;
[randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_multi_sim_theory_multi(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows,numKymos,numTheories);


% [randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_multi_sim_theory_single(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows,numKymos);
% 
% folder with kymos
dataFold = sets.outFold;
% folder with theories

% now run comparison
numFrames =1;
nmBp = 0.3;
comparisonMethod = 'mpAll';
% comparisonMethod = 'mp';
% comparisonMethod = 'spearman';
% comparisonMethod = 'dtw';
[rezMax,bestBarStretch,bestLength,rezMaxAll] = run_parallel_str('hca_p1.txt', theoriesTxt, dataFold, numFrames, nmBp,comparisonMethod);

thrI = 1;
barI = 2;
strI = 4;
posVals = rezMaxAll{thrI}{barI}{strI}.pos;

figure,plot(posVals)
figure,plot(rezMaxAll{thrI}{barI}{strI}.maxcoef)
% we can compute how much is the data stretched. This is a good alternative
% to getting the correct stretching factor, as it is based on all
% subfragments (rather than the best subfragment)
params = polyfit(1:length(posVals),posVals', 1);
params(1)

figure,hold on
for i=1:length( rezMaxAll{thrI}{barI})
    
    posVals = rezMaxAll{thrI}{barI}{i}.pos;
    params = polyfit(1:length(posVals),posVals', 1);
    if abs(params(1)) < 0.1
            plot( rezMaxAll{thrI}{barI}{i}.pos,'x')
    end
% figure,plot(posVals)
% figure,plot(rezMaxAll{thrI}{barI}{strI}.maxcoef)
% % we can compute how much is the data stretched. This is a good alternative
% to getting the correct stretching factor, as it is based on all
% subfragments (rather than the best subfragment)



end