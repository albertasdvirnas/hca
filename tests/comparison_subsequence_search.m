
sets = ini2struct( 'sets.txt' );



seq2 = [];
for idx=1:10
    seq2{idx} = imgaussfilt(normrnd(0,1,1,500000),3);
end
mkdir('resultData');

theoryStruct = save_rand(seq2,'resultData/');

seq1 = [];
idx=1;
% generate one 
seq1{idx} = imgaussfilt(normrnd(0,1,1,500),3);
barcodeGen{idx}.rawBarcode = seq1{idx};
barcodeGen{idx}.rawBitmask = ones(1,length(seq1{idx}));
barcodeGen{idx}.rawBitmask(1:5) = 0;
barcodeGen{idx}.rawBitmask(end-7:end) = 0;
idx2=1;

fileID = fopen(theoryStruct{idx2}.filename,'r');
formatSpec = '%f';
theorBar = transpose(fscanf(fileID,formatSpec));
fclose(fileID);

    
import SignalRegistration.unmasked_pcc_corr;
tic
xcorrsOld = unmasked_pcc_corr(barcodeGen{idx}.rawBarcode, theorBar, barcodeGen{idx}.rawBitmask);
toc

% TODO: for large k, we don't get the same values


tic
xcorrs2 = unmasked_MASS_PCC(theorBar, barcodeGen{idx}.rawBarcode,barcodeGen{idx}.rawBitmask, 2^12);

% xcorrs2 = MASS_PCC([theorBar theorBar(1:length(barcodeGen{idx}.rawBarcode)-1)], barcodeGen{idx}.rawBarcode, 2^9);
toc

% theorBar
% save this plot to show that the results are identical up to numerical
% precision.
figure,plot((xcorrsOld(:)-xcorrs2(:)))
% [rezMax{idx}.maxcoef,rezMax{idx}.pos,rezMax{idx}.or] = get_best_parameters(xcorrs, 3 );

xcorrs = MASS_PCC([theorBar theorBar(1:length(barcodeGen{idx}.rawBarcode)-1)], barcodeGen{idx}.rawBarcode, 2^9);

 sum((xcorrs2(:)-xcorrsOld(:)).^2)
vals = [];
kk=9:22;
for k=kk
    xcorrs2 = MASS_PCC([theorBar theorBar(1:length(barcodeGen{idx}.rawBarcode)-1)], barcodeGen{idx}.rawBarcode, 2^k);
    vals = [vals sum((xcorrs(:)-xcorrs2(:)).^2)];
end

figure
hold on
plot(kk,vals)
hold on
plot(kk(end), sum((xcorrs(:)-xcorrsOld(:)).^2),'redx') 


%% how would it work on bacterial stuff:

xx = normrnd(0,1,1,125823789);
y = normrnd(0,1,1,500);

% we could save just maximum
tic
xcorrs2 = MASS_PCC(xx,y,2^(2+nextpow2(length(y))));
toc