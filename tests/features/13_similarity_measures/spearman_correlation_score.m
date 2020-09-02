% test rand kymo scores for each row when aligning to theory

sets =ini2struct('sp_sim_sets_1.txt');
setstxt = 'sp_settings.txt';

% % chech if the sequences have to be reproducible
rng(1,'twister');
    
numRows = 5;
 sets.outFold = 'test/';
[randKymos,names,posCuts,lengths,theoriesTxt,ver] = generate_sim_theory_single(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows);
% 
dataFold = sets.outFold;

% settings
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings(setstxt);
sets.whichtokeep = 1;
    
listing = dir(fullfile(dataFold,'*.tif'));

% put the names in a txt file, could already have a txt file.
fd =fopen('tifs.txt','w');
for i=1:length(listing)
    fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
end
fclose(fd);
sets.kymosets.kymoFile = 'tifs.txt';
        
sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


import CBT.Hca.Settings.get_user_settings;
sets = get_user_settings(sets);
try mkdir(sets.output.matDirpath);catch; end
% add kymographs
import CBT.Hca.Import.add_kymographs_fun;
[kymoStructs] = add_kymographs_fun(sets);

%  put the kymographs into the structure
import CBT.Hca.Core.edit_kymographs_fun;
kymoStructs = edit_kymographs_fun(kymoStructs,0);

import generate.kymo_to_multi_bar;
kymoStructs = kymo_to_multi_bar(kymoStructs);
% now convert this to many kymo of single frame

% align kymos - should not take any time if it's just single frame
import CBT.Hca.Core.align_kymos;
[kymoStructsAligned] = align_kymos(sets,kymoStructs);
      
   % generate barcodes
import CBT.Hca.Core.gen_barcodes;
barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructsAligned, sets);
 
trueStart = -ver+50;
startFound = cellfun(@(x) x.lE,barcodeGen);

% %% !!!
% % this shows if we estimate the start pixel correctly (based on PSF)
% figure,plot(startFound);hold on;plot(trueStart); legend({'Found','True'})
% figure,plot(startFound-trueStart)

sets.theories = theoriesTxt;
 sets.theory.nmbp = 0.3;
% get user theory
import CBT.Hca.Settings.get_user_theory;
[theoryStruct, sets] = get_user_theory(sets);

%%
% len1 = 500;
% len2 = 500;
% psf = 3;
% % define query and data 
% A = imgaussfilt(rand(len1,1),psf);
% B = imgaussfilt(rand(len2,1),psf);
% 
% W1 = ones(len1,1);
% 
% % compute ranks of A and B
% Ar = compute_ranks(A);
% Br = compute_ranks(B);
% pcc = @(x,y) zscore(x,1)*zscore(y',1)/length(x);
% pcc(Ar',Br')
% 
% [RHO,PVAL] = corr(A,B,'Type','Spearman');
% [RHO,PVAL] = corr(Ar,Br,'Type','Pearson');
% 
% 
% % [Ar,ArIdx] = sort(A);
% 
% import SignalRegistration.unmasked_spearman_corr;
%  [maxSpearman,posMax,orMax,scoreVec] = unmasked_spearman_corr(A,B,W1);
 
 % A = [1 4 6 3 4 6 7 8]; 
% b = [34 56 34 56 79 23 48 28];
% [RHO,PVAL] = corr(a',b','Type','Spearman');
% query = 1:10;
% data = 1:10;
% data(4) = 5;


%% The method should be run through on_compare function

% note: this might take a lot of space for big B, so might not be
% efficient..
%pieces with the length of k and with overlaps of m-1 
%samples and make a matrix s from the data vector. 
% s=buffer(B,length(A),length(A)-1); 
            
            
            
sets.w = 200;
sets.comparisonMethod = 'spearman';
tic
import CBT.Hca.Core.Comparison.on_compare;
[rezMax,bestBarStretch,bestLength] = on_compare(barcodeGen,theoryStruct{1},sets.comparisonMethod,0.98,sets.w);

% [rezMax,bestBarStretch,bestLength] = on_compare(barcodeGen,theoryStruct{1},sets.comparisonMethod,sets.theory.stretchFactors,sets.w);
toc


   %% other methods
sets.w = 200;
sets.comparisonMethod = 'mp';
tic
import CBT.Hca.Core.Comparison.on_compare;
[rezMax,bestBarStretch,bestLength] = on_compare(barcodeGen,theoryStruct{1},sets.comparisonMethod,1,sets.w);
toc