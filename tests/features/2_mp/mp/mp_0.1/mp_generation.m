% We want to use mp in
% [ comparisonStructure ] = on_compare_theory_to_exp( barcodeGen,theoryStruct, sets)
addpath(genpath(pwd));

len1 = 500; % length query
len2 = 10000; % length data
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3); % simulated query
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3); % simulated data

bit1 = ones(len1,1); % bitmask, assumed to be ones

r = 300;
k= 2^14;
maskWindow = 50;


import mp.mp_masked_profile_stomp_dna;
comparisonFun = @(x,y,z,w) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))));

tic
[mp, mpI] = comparisonFun(bar1,bar2,bit1,r);
toc

[ maxcoef,pos,or ] = get_best_parameters_mp( mp,mpI, 3, len1,len2, 1, maskWindow);

