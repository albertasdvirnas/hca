% test for pcc for BAC subfragments


% timestamp for the results
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% import settings
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings('bac_settings.txt');
sets.resultsDir = 'out/';
mkdir(sets.resultsDir);

% preparation: data
% J21
data = '/home/albyback/workData/rawData/HumanDNAProject/BACs/Kymographs/J21/';
theory = '/home/albyback/workData/rawData/HumanDNAProject/BACs/Sequences/J21.fasta';
sets.theoryGen.meanBpExt_nm = 0.223;

% % P18
% data = '/media/albyback/My Passport/workData/HC/Human DNA Project/BACs/Kymographs/P18/';
% theory = '/media/albyback/My Passport/workData/HC/Human DNA Project/BACs/Sequences/P18 (with additional sequence one side).fa';
% sets.theoryGen.meanBpExt_nm = 0.225;

% %H8
% data = '/media/albyback/My Passport/workData/HC/Human DNA Project/BACs/Kymographs/H8/';
% theory = '/media/albyback/My Passport/workData/HC/Human DNA Project/BACs/Sequences/H8 (with surrounding sequence).fa';
% sets.theoryGen.meanBpExt_nm = 0.235;

sets.timeFramesNr = 10;
sets.alignMethod = 3;
% gen barcodes
barcodeGen = gen_bac(data,sets);
 
% gen theory
[theoryGen] = gen_theory(theory,sets);

%% now run MP:

 sets.k = 2^9;
sets.c = 200;
sets.stretch = 0.95:0.01:1.05;
%
maxCcof = zeros(1,length( barcodeGen));
comparisonStruct = cell(1,length( barcodeGen));
for i=1:length( barcodeGen)
%     i
    barcode = barcodeGen{i}.rawBarcode( barcodeGen{i}.rawBitmask);

    
    import mp.mp_stretch;
    [comparisonStruct{i}] = mp_stretch(barcode,theoryGen.theoryBarcodes{1}, sets.c,sets.k,sets.stretch);
    maxcc = cellfun(@(x) max([0; x.mp]), comparisonStruct{i});
    maxCcof(i) = max(maxcc);
end

maxCcof(maxCcof==0)= nan;
f =figure,plot(maxCcof,'x'); hold on, plot(repmat(nanmean(maxCcof),1,length(maxCcof)))
legend({'Individual max scores','average'},'location','best')

mean(maxCcof)


[a,b,c]=fileparts(theory)
title(b)
saveas(f,fullfile(sets.resultsDir,strcat(b,'maxfragmentscore.eps')),'epsc')
saveas(f,fullfile(sets.resultsDir,strcat(b,'maxfragmentscore.png')))

% 
% figure,plot(comparisonStruct{i}{3}.mpI)
%     load(fullfile(listing(i).folder,listing(i).name));
% 
% bar1 =    imresize(barcode', [1,round(length(barcode)*0.95)]);
% bar2 = theoryGen.theoryBarcodes{1};
% % 

%%
i=23
maxcc = cellfun(@(x) max(x.mp), comparisonStruct{i});
[c,d] = max(maxcc);
sets.stretch(d)
idy = d;

barcode = barcodeGen{i}.rawBarcode( barcodeGen{i}.rawBitmask);

bar1 =    imresize(barcode, [1,round(length(barcode)*sets.stretch(d))]);
bar2 = theoryGen.theoryBarcodes{1};
% 

[a,b] = max(comparisonStruct{i}{idy}.mp)
b1 = zscore(bar1(b:b+sets.c-1));
st = comparisonStruct{i}{idy}.mpI(b);
sto = comparisonStruct{i}{idy}.mpI(b)+sets.c-1;
if st > length(bar2)
    st = mod(st-1,length(bar2))+1;
    sto = st + sets.c-1;
    b2 = zscore(fliplr(bar2(st:sto)));
else
    b2 = zscore(bar2(st:sto));
end


f=figure,plot(b1)
hold on
plot(b2)
[~,b,c]=fileparts(theory)
% plot(-49:length(b1)+50,zscore(bar2(st-50:sto+50))-1,'black')
legend({'experiment',strcat('Theory_',b)})
title(strcat(['Score=' num2str(a)  ]))

saveas(f,fullfile(sets.resultsDir,strcat(b,'bestcomparison.png')))
