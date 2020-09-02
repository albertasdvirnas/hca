% long over-due feature - combining SV to HCA..

% I guess if we have barcodeGen structure, theoryStruct and corresponding
% settings, then we could just run the HMM module. How best to integrate
% HMM module with this?

% call the class hmm_search for now.

% as input, we have barcodeGen, theoryStruct and settings, at first we
% simulate all this here.


sets =ini2struct('sp_sim_sets_1.txt');
setstxt = 'sp_settings.txt';
% % chech if the sequences have to be reproducible
rng(1,'twister');
    
% could put this to hca_search too
numRows = 10;
 sets.outFold = 'test/';
 
obj =  hca_search.generate_sim('sp_sim_sets_1.txt',10,'test/');
%  
% [randKymos,names,posCuts,lengths,theoriesTxt,ver,tifsTxt] = generate_sim_theory_single(sets.lengthN,sets.kernelsigma,sets.outFold,sets.lambda,sets.pccScore,sets.strF,numRows);
% % 
% dataFold = sets.outFold;
% sets.kymosets.kymoFile = tifsTxt;

%% now hca stuff
% setstxt = 'sp_settings.txt';

% hca structure
obj = hca_search( 'sp_settings.txt',obj.theoriesTxt,obj.tifsTxt,'out/');

% this can be done also for chromosome fragments which were already matched
% to best theory
fig1=figure;
idx = 5;
splitStr = strsplit(obj.barcodeGen{idx}.name,'_');
comparisonStruct{idx}.pos = str2num(splitStr{2});
   
comparisonStruct{idx}.idx = 1;
comparisonStruct{idx}.or = 1;
comparisonStruct{idx}.bestBarStretch =  1/str2num(splitStr{4});
import CBT.Hca.UI.Helper.plot_best_bar;
maxcoff = zeros(10,3);
maxcoff(idx,1) = 1;
plot_best_bar(fig1,obj.barcodeGen,[],comparisonStruct, obj.theoryStruct, maxcoff,100);


import Plot.alignment_class;
sets.fold = 'out';
alignment_class(obj.barcodeGen{1},obj.theoryStruct{1},sets,'bla')



% now just need to run SV
% 
% % settings
% import CBT.Hca.Import.import_hca_settings;
% [sets] = import_hca_settings(setstxt);
% 
% 
% sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% sets.kymosets.kymoFile = tifsTxt;
% 
% 
% import CBT.Hca.Settings.get_user_settings;
% sets = get_user_settings(sets);
% sets.output.matDirpath = 'out/';
% sets.whichtokeep = 1;
% try mkdir(sets.output.matDirpath);catch; end
% 
% % add kymographs
% import CBT.Hca.Import.add_kymographs_fun;
% [kymoStructs] = add_kymographs_fun(sets);
% 
% %  put the kymographs into the structure
% import CBT.Hca.Core.edit_kymographs_fun;
% kymoStructs = edit_kymographs_fun(kymoStructs,0);
% 
% import generate.kymo_to_multi_bar;
% kymoStructs = kymo_to_multi_bar(kymoStructs);
% % now convert this to many kymo of single frame
% 
% % align kymos - should not take any time if it's just single frame
% import CBT.Hca.Core.align_kymos;
% [kymoStructsAligned] = align_kymos(sets,kymoStructs);
%       
%    % generate barcodes
% import CBT.Hca.Core.gen_barcodes;
% barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructsAligned, sets);
% %  
% % trueStart = -ver+50;
% % startFound = cellfun(@(x) x.lE,barcodeGen);
% 
% % %% !!!
% % % this shows if we estimate the start pixel correctly (based on PSF)
% % figure,plot(startFound);hold on;plot(trueStart); legend({'Found','True'})
% % figure,plot(startFound-trueStart)
% 
% sets.theories = theoriesTxt;
%  sets.theory.nmbp = 0.3;
% % get user theory
% import CBT.Hca.Settings.get_user_theory;
% [theoryStruct, sets] = get_user_theory(sets);



