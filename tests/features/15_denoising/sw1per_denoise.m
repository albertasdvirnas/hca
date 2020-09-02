
sets =ini2struct('sp_sim_sets_1.txt');
setstxt = 'sp_settings.txt';

% % chech if the sequences have to be reproducible
rng(1,'twister');
    
numRows = 1;
numKymos = 20;
numTheories = 1;
sets.lengthN = repmat(10000,1,numTheories);
sets.outFold = 'test/';
sets.strF = 2;
sets.pccScore = 0.9;
[randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_multi_sim_theory_multi(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows,numKymos,numTheories);

% [randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_multi_sim_theory_single(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows,numKymos);
% 
% folder with kymos
dataFold = sets.outFold;
% folder with theories

% now run comparison
numFrames =1;
nmBp = 0.3;
comparisonMethod = 'mass_pcc';
comparisonMethod = 'mp';
comparisonMethod = 'spearman';
[hcaStruct] = run_parallel('hca_parallel_settings.txt', theoriesTxt, dataFold, numFrames, nmBp,comparisonMethod);
% [hcaStruct] = run_parallel('hca_parallel_settings.txt', theoriesTxt, dataFold, numFrames,hcaStruct)

% randKymos
bar = randKymos{1};
bar = movmean(bar,2)
splineF = @(x) interp1(1:length(bar),bar,x);

figure,plot(splineF(1.2:1:length(bar)-1))