% Consider some kymograph, i.e. P18
% import settings
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings('bac_settings.txt');
sets.resultsDir = 'out/';
mkdir(sets.resultsDir);

data = dir("/media/albyback/My Passport/DATA/Shared - Kymographs for Albertas testing/Ecoli_BL21/*/*.tif");
kymos = fullfile(data(:).folder, data(:).name);

sets.names = arrayfun(@(x) x.name, data,'UniformOutput',false);
sets.folds = arrayfun(@(x) x.folder, data,'UniformOutput',false);

idx = 2;

A = imread(fullfile(sets.folds{idx},sets.names{idx}));

import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
[ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(A), sets.edgeDetectionSettings);
% 
% xx = 15;
% figure,plot(A(xx,leftEdgeIdxs(xx): rightEdgeIdxs(xx)));
% hold on
% plot(A(xx+1,leftEdgeIdxs(xx): rightEdgeIdxs(xx)));

% now, we want to find a best matching subfragment between each row, we can
% use MP


tic
sets.k = 2^9;
sets.c = 300;
% sets.stretch = 0.95:0.01:1.05;
sets.stretch = 1;
%
maxCcof = zeros(1,length( barcodeGen));
maxStretch = zeros(1,length( barcodeGen));

comparisonStruct = cell(1,length( barcodeGen));
posI = zeros(1,size(A,1)-1);
posJ = zeros(1,size(A,1)-1);

fullkymo = nan(2*size(A,1),size(A,2));
for i=1:size(A,1)-1
%     i
    barcode1 = double(A(i,:));
    barcode2 = double(A(i+1,:));
    
    % take just the barcode indices
    barcode1 = barcode1(leftEdgeIdxs(i):rightEdgeIdxs(i));
    % take just the barcode indexes
    barcode2 = barcode2(leftEdgeIdxs(i+1):rightEdgeIdxs(i+1));
    
    
    % todo: this case mpI should not include circular stuff or flipped
    % case, so it is just the regular MP!
    
    % compute the full mp profile
    import mp.mp_stretch_regular;
    [comparisonStruct{i}] = mp_stretch_regular(barcode1,barcode2, sets.c,sets.k,sets.stretch);
    maxcc = cellfun(@(x) max([0; x.mp]), comparisonStruct{i});
    [maxCcof(i),maxStretch(i)] = max(maxcc);
    
    % get positions on a and b
    [a,positionBb1] = max(comparisonStruct{i}{maxStretch(i)}.mp);
    positionBb2 = comparisonStruct{i}{maxStretch(i)}.mpI(positionBb1);
    
    % get these back to original positions
    posI(i) = positionBb1+leftEdgeIdxs(i)-1;
    posJ(i) = positionBb2+leftEdgeIdxs(i+1)-1;

    barcode1 =   imresize(barcode1, [1,round(length(barcode1)*sets.stretch(maxStretch(i)))]);
    fullkymo(2*i-1,posI(i):posI(i)+sets.c-1) = barcode1(positionBb1:positionBb1+sets.c-1);
    fullkymo(2*i,posJ(i):posJ(i)+sets.c-1) = barcode2(positionBb2:positionBb2+sets.c-1);

end
toc


%% kymo alignment to first row



tic
sets.k = 2^9;
sets.c = 200;
sets.stretch = 0.9:0.01:1.1;
% sets.stretch = 1;
%
maxCcof = nan(1,length( barcodeGen));
maxStretch = zeros(1,length( barcodeGen));

comparisonStruct = cell(1,length( barcodeGen));
posI = zeros(1,size(A,1)-1);
posJ = zeros(1,size(A,1)-1);

fullkymo = nan(size(A,1),size(A,2));

barcode1 = double(A(1,:));
% take just the barcode indices
barcode1 = barcode1(leftEdgeIdxs(1):rightEdgeIdxs(1));

    
for i=2:size(A,1)-1
%     i
    barcode2 = double(A(i+1,:));
       % take just the barcode indexes
    barcode2 = barcode2(leftEdgeIdxs(i+1):rightEdgeIdxs(i+1));
    
    
    % todo: this case mpI should not include circular stuff or flipped
    % case, so it is just the regular MP!
    
    % compute the full mp profile
    import mp.mp_stretch_regular;
    [comparisonStruct{i}] = mp_stretch_regular(barcode1,barcode2, sets.c,sets.k,sets.stretch);
    maxcc = cellfun(@(x) max([0; x.mp]), comparisonStruct{i});
    [maxCcof(i),maxStretch(i)] = max(maxcc);
    
    % get positions on a and b
    [a,positionBb1] = max(comparisonStruct{i}{maxStretch(i)}.mp);
    positionBb2 = comparisonStruct{i}{maxStretch(i)}.mpI(positionBb1);
    
    % get these back to original positions
    posI(i) = positionBb1+leftEdgeIdxs(i)-1;
    posJ(i) = positionBb2+leftEdgeIdxs(i+1)-1;

    barcode1 =   imresize(barcode1, [1,round(length(barcode1)*sets.stretch(maxStretch(i)))]);
%     fullkymo(i,posI(i):posI(i)+sets.c-1) = barcode1(positionBb1:positionBb1+sets.c-1);
    fullkymo(i,posJ(i):posJ(i)+sets.c-1) = barcode2(positionBb2:positionBb2+sets.c-1);

end
toc


% 
% isLinearTF =1;
% mask = 10;
% import CBT.Hca.UI.Helper.get_best_parameters_mp;
% [ maxcoef,pos,or ] = get_best_parameters_mp( mp,mpI, 10, len1,len2,  isLinearTF, mask);
% 

% sets.edgeDetectionSettings.method = 'Error function';
% import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
% [ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(A), sets.edgeDetectionSettings);


% first find edges and bitmask, so that we wouldn't go out


