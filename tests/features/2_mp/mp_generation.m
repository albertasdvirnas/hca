% We want to use mp in
% [ comparisonStructure ] = on_compare_theory_to_exp( barcodeGen,theoryStruct, sets)

len1 = 500;
len2 = 1000;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
% bar2 = [imgaussfilt(normrnd(0,1,len2,1),2.3); bar1; imgaussfilt(normrnd(0,1,len2,1),2.3)];

bit1 = ones(len1,1);

r = 300;
k= 2^14;


import mp.mp_masked_profile_stomp_dna;
comparisonFun = @(x,y,z,w) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))));

tic
[mp, mpI] = comparisonFun(bar1,bar2,bit1,r);
toc

import CBT.Hca.UI.Helper.get_best_parameters_mp;
[ maxcoef,pos,or ] = get_best_parameters_mp( mp,mpI, 3, len1,len2, 1, 50,r);

% check that this actually gives the correct position
figure,plot(bar2)
hold on
plot(pos(1):pos(1)+length(bar1)-1,bar1)
frag1 = flipud(bar1);
frag2 = bar2(pos(1):pos(1)+length(bar1)-1);
pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
pcc(frag1,frag2)

figure,plot(fliplr(frag1))
hold on
plot(frag2)
% 
% import CBT.Hca.UI.Helper.get_best_parameters;
% [maxcoef,pos,orient] = get_best_parameters(xcorrs, 3,length(barC),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);



comparisonFun = @(x,y,z) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));
tic
xcorrs = comparisonFun(bar1',bar2',bit1');
toc

% 
% import mp.mp_profile_stomp_dna;
% comparisonFun = @(x,y,w) mp_profile_stomp_dna(y,x,2^(4+nextpow2(length(x))));

%% mp gives the subframent on barcode 