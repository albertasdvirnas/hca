% First simulate some kymoStruct
load('cbc_test.mat');

sets.length = 300;
sets.kernelsigma = 2.3;
sets.rand.randmethod = 'imgaussfilt';
sets.pccScore = 0.8;
 bar1 = normrnd(0,1,1,sets.length );
 
bar1 =  1+imgaussfilt(bar1,sets.kernelsigma);

import Rand.add_noise_to_barcode;

numFrames = 20;
kymo = zeros(numFrames,500);

strF =  normrnd(0,0.01,1,numFrames );

stPos = randi([100 150],numFrames);

for i=1:numFrames
    [bar1N] = add_noise_to_barcode(bar1, sets);
    kymo(i,:) = exprnd(1/10,1,500);
    newLen = round(sets.length*(1+strF(i)));
    kymo(i,stPos(i):stPos(i)+newLen-1) = imresize(bar1N,[1,newLen]);
end

kymoStructs{1}.unalignedKymo = kymo;
figure,imagesc(kymo)


import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
[ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(kymo, sets.edgeDetectionSettings);

idx = 1;
kymoSingleFrameStruct = [];
sz = size(kymoStructs{idx}.unalignedKymo);
kymos =mat2cell(kymoStructs{idx}.unalignedKymo,ones(1,sz(1)),sz(2));

if sets.timeFramesNr == 0
    minV = inf;
else
    minV = sets.timeFramesNr;
end

i=1;
% for i=1:min(minV,size(kymoStructs{idx}.unalignedKymo,1))
    kymoSingleFrameStruct{i}.unalignedKymo = kymo;

    kymoSingleFrameStruct{i}.alignedKymo = double(kymos{i});
    kymoSingleFrameStruct{i}.name = kymoStructs{idx}.name;
%     kymoSingleFrameStruct{i}.leftEdgeIdxs = leftEdgeIdxs(i);
%     kymoSingleFrameStruct{i}.rightEdgeIdxs = rightEdgeIdxs(i);
% % end
    % 
    
    
import CBT.Hca.Core.align_kymos_consensus;
[ kymoStructs,barcodeGenMol,consensusStructs ] = align_kymos_consensus( sets, kymoSingleFrameStruct );




% align kymos
import CBT.Hca.Core.align_kymos;
[kymoStructs] = align_kymos(sets,kymoStructs);
