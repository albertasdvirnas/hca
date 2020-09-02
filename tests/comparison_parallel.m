% Testing speed of parallel comparison
testTheories = 'paralleltesting.txt';
import CBT.Hca.Import.import_names;
[names] = import_names(testTheories);

% theories = load(names{1});

% Some random experiment
testExperiment = 'parallelexperiment.txt';
[namesExperiment] = import_names(testExperiment);

% theoryStruct = save_rand(theories.hcaSessionStruct.theoryGen.theoryBarcodes,'resultData2/');

theoryStruct = [];
listing = dir('resultData2/*.txt');
names = {listing(:).name};
for i=1:min(10000,length(names))
    theoryStruct{i}.filename = names{i};
end

%% 
sets = ini2struct( 'sets.txt' );


seq1 = imgaussfilt(normrnd(0,1,1,500),3);

seq2 = [];
for idx=1:10
    seq2{idx} = imgaussfilt(normrnd(0,1,1,50000),3);
end
mkdir('resultData');

% theoryStruct = save_rand(seq2,'resultData/');

% generate one 
barcodeGen{1}.rawBarcode = seq1;
barcodeGen{1}.rawBitmask = ones(1,length(seq1));

comparisonStruct1 = cell(1,length(theoryStruct));
comparisonStruct2 = cell(1,length(theoryStruct));
comparisonStruct3 = cell(1,length(theoryStruct));


sets.comparisonMethod = 'unmasked_pcc_corr';

% ideal: barcodeGen is available (sent to all the workers, and each worker
% gets part of theoryStruct (ideally so that total length would be roughly
% similar..)
tic
% unfiltered comparison
parfor barNr = 1:length(theoryStruct)
%     disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));

%     import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
    comparisonStruct1{barNr} = CBT.Hca.Core.Comparison.on_compare_theory_to_exp(barcodeGen,theoryStruct{barNr}, sets);
end
toc


sets.comparisonMethod = 'unmasked_pcc_corr';

tic
% unfiltered comparison
for barNr = 1:length(theoryStruct)
%     disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));

%     import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
    comparisonStruct2{barNr} = CBT.Hca.Core.Comparison.on_compare_theory_to_exp(barcodeGen,theoryStruct{barNr}, sets);
end
toc


sets.comparisonMethod = 'mass_pcc';

tic
sets.theory.stretchFactors = 1;
% unfiltered comparison
for barNr = 1:length(theoryStruct)
%     disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));

%     import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
    comparisonStruct3{barNr} = CBT.Hca.Core.Comparison.on_compare_theory_to_exp(barcodeGen,theoryStruct{barNr}, sets);
end
toc
%     
%     
% tic
% import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
% [ comparisonStructure ] = on_compare_theory_to_exp( barcodeGen,theoryStruct, sets);
% toc