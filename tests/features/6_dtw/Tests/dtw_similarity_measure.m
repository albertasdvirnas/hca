% this code to test DTW similarity measure with simulated data
% first generate some data

sets =ini2struct('sp_sim_sets_1.txt');
setstxt = 'sp_settings.txt';

% % chech if the sequences have to be reproducible
rng(1,'twister');
    
% one kymo and one theory
numRows = 3;
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
% dataFold = sets.outFold;
% folder with theories

% now run comparison
numFrames =1;
nmBp = 0.3;
comparisonMethod = 'dtw';
% comparisonMethod = 'mp';
% comparisonMethod = 'spearman';
% comparisonMethod = 'dtw';

% 
dataFold = sets.outFold;

% settings
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings(setstxt);
sets.whichtokeep = 1:2;
    
listing = dir(fullfile(dataFold,'*.tif'));

% put the names in a txt file, could already have a txt file.
fd =fopen('tifs.txt','w');
for i=1:length(listing)
    fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
end
fclose(fd);
sets.kymosets.kymoFile = 'tifs.txt';
        
sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sets.output.matDirpath = 'output';

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
 
trueStart = cell2mat(cellfun(@(x) arrayfun(@(y) -y+50,x),ver,'UniformOutput',false));
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

            
sets.w = 200;
sets.comparisonMethod = 'dtw';
% tic
import CBT.Hca.Core.Comparison.on_compare;
[rezMax,bestBarStretch,bestLength] = on_compare(barcodeGen,theoryStruct{1},sets.comparisonMethod,0.98,sets.w);
% 


% [rezMax,bestBarStretch,bestLength,rezMaxAll] = run_parallel_str('hca_p1.txt', theoriesTxt, dataFold, numFrames, nmBp,comparisonMethod);
