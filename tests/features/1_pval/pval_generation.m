% First generate some random series

% want to generate p-values for bacterial chromosomes

load('/home/albyback/data/Ref to Albertas/Ref_Theory_ncbi_chromosomes_190717_130nmPixel_0.34nmPERbp_HCAv3.8.3.mat');

% lengths of theory barcodes
lengths = cellfun(@(x) length(x), hcaSessionStruct.theoryGen.theoryBarcodes);

totalLength = sum(lengths); % 25 times more than HCA.

sets = ini2struct('pval_sets.txt');
sets.pccScore = 0.95;
sets.theory.isLinearTF = 1;
sets.length=200;
sets.numRnd = 1000;
sets.length2 = 50;
iiVals = 0:36;
lenQ = 200;
sets.len2 = sets.length;
sets.fullPath = 'pval.txt';
sets.islinearQ = 1; % linear query/not true for consensus
sets.islinearD = 1; % linear data/ humans have linear, but bacteria circular
sets.lenMax = 110;
sets.lenMin = 100;
sets.lengths = lengths(1:1000);


import CBT.Hca.Core.Pvalue.precompute_pvalue_files_multi_theory;
dataCC = precompute_pvalue_files_multi_theory('pval.txt', sets );

load('Development/hca/tests/features/1_pval/lengths.mat');
[dataCC,t] = CBT.Hca.Core.Pvalue.precompute_pvalue_files_multi_theory('pval_sets.txt',lengths(1:5));
