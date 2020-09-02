% first generate some data, so that kymo's and theory could be used as
% input

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
% comparisonMethod = 'mp';
% comparisonMethod = 'spearman';
% comparisonMethod = 'dtw';
[hcaStruct] = run_parallel('hca_parallel_settings2.txt', theoriesTxt, dataFold, numFrames, nmBp,comparisonMethod);
% [hcaStruct] = run_parallel('hca_parallel_settings.txt', theoriesTxt, dataFold, numFrames,hcaStruct)

% [hcaStruct] = HCA_run_parallel('hca_parallel_settings.txt', theoriesTxt, dataFold, 1);

% t = HCA_theory_parallel('/home/albyback/workData/rawData/HumanDNAProject/BACs/Sequences/*.fa',0.225);
