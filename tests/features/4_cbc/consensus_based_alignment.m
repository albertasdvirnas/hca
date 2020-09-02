load('cbc_test.mat');

% we can use consensus code to align a kymograph - this does stretch each
% row to the same length (which should be ok since we detect the edges),
% the rows which are completely off will have a very low CC score

% for kymographs, when averaging use original intensities
sets.consensus.barcodeNormalization= 'original';
sets.timeFramesNr = 10;

% take some index / or loop through all, might be slower
for idx = 1:size(kymoStructs,2);

    % we use kymoStructs to split each kymo to it's separate kymostruct, then
    % we can generate barcodes. We save each row of unaligned kymograph as
    % aligned kymo. Need to compute edges for the unaligned kymo (instead of
    % aligned)
    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(kymoStructs{idx}.unalignedKymo), sets.edgeDetectionSettings);

    kymoSingleFrameStruct = [];
    sz = size(kymoStructs{idx}.unalignedKymo);
    kymos =mat2cell(kymoStructs{idx}.unalignedKymo,ones(1,sz(1)),sz(2));

    if sets.timeFramesNr == 0
        minV = inf;
    else
        minV = sets.timeFramesNr;
    end

    for i=1:min(minV,size(kymoStructs{idx}.unalignedKymo,1))
        kymoSingleFrameStruct{i}.alignedKymo = double(kymos{i});
        kymoSingleFrameStruct{i}.name = kymoStructs{idx}.name;
        kymoSingleFrameStruct{i}.leftEdgeIdxs = leftEdgeIdxs(i);
        kymoSingleFrameStruct{i}.rightEdgeIdxs = rightEdgeIdxs(i);
    end
    % 

    % now we can generate barcodes
    import CBT.Hca.Core.gen_barcodes;
    barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoSingleFrameStruct, sets);

    % barcode stretching factors. Note that edges are calculated only
    % approximately, but some tendencies in the stretchings can be seen here.
    % i.e stretching not completely random in between rows 
    % f = figure,plot(cellfun(@(x) x.stretchFactor,barcodeGen));
    % xlabel("Timeframe nr.");
    % ylabel("Stretch factor");
    % title("Stretch factors plot");
    % legend({kymoSingleFrameStruct{1}.name},'Interpreter','latex','location','best');
    sets.output = 'out';
    % saveas(f,fullfile('out','stretching_factors.eps'))




    % we want to use generate consensus, theerefore we need to plug in
    % barcodeGen structure
    % % generate consensus
    % import CBT.Hca.Core.gen_consensus;
    % consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);

    % note that in this case barcodes are linear, and can only have a few pixel
    % shifts in between them, so instead of full FFT, it would suffice to run
    % PCC in the neighborhood (or a linear PCC), which would require us adding
    %     linear=1;
    %     if linear
    %         xcorrs(1,(length(longVec)-sum(w1)+1):end)=0;
    %         xcorrs(2,(length(longVec)-sum(w1)+1):end)=0;
    %     end
    % to the cc calculation
    import CBT.Hca.Core.gen_consensus;
    consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);


    import CBT.Hca.UI.Helper.select_all_consensus;
    [consensusStruct,sets] = select_all_consensus(consensusStructs,sets);

    consensusStructs.consensusStruct = consensusStruct;

    hcaStruct = [];
    hcaStruct.barcodeGenC = barcodeGen;
    hcaStruct.sets = sets;
    hcaStruct.consensusStruct = consensusStruct;
    hcaStruct.consensusStructs = consensusStructs;
    hcaStruct.kymoStructs = kymoStructs;

    % import CBT.Hca.Export.export_cbc_compatible_consensus;
    % export_cbc_compatible_consensus(consensusStructs,barcodeGen,kymoStructs,sets);

    % if we want to plot
%     import CBT.Hca.Export.plot_consensus_concentrically;
%     plot_consensus_concentrically(consensusStructs,barcodeGen)

    kymoStructs{idx}.consensusStructs = consensusStructs;
     kymoStructs{idx}.barcodeGen = barcodeGen;
    kymoStructs{idx}.alignedKymoCBC = consensusStructs.treeStruct.barcodes{end};
    kymoStructs{idx}.barcode = consensusStructs.consensusStruct{end}.rawBarcode;
    kymoStructs{idx}.bitmask = consensusStructs.consensusStruct{end}.rawBitmask;
    barcodeGenMol{idx}.bgMeanApprox = mean(cellfun(@(x) x.bgMeanApprox,barcodeGen));
    barcodeGenMol{idx}.bgStdApprox =  mean(cellfun(@(x) x.bgStdApprox,barcodeGen));
    barcodeGenMol{idx}.rawBarcode = consensusStructs.consensusStruct{end}.rawBarcode;
    barcodeGenMol{idx}.rawBitmask = consensusStructs.consensusStruct{end}.rawBitmask;
%     barcodeGenMol{idx}.lE = 

end

numBar = length(barcodeGenMol);
sets.consensus.barcodeNormalization= 'original';

%% Prestretch barcodes to the same length
% could also try not to stretch, or use the 5% stretch?
% convert to common length, if chosen
if  sets.genConsensus == 1
    allLengths = cellfun(@(x) length(x.rawBarcode),barcodeGenMol);
    commonLength = ceil(mean(allLengths));
    stretchings = commonLength./allLengths;
    strMin = min(stretchings);
    strMax =  max(stretchings);
    disp(strcat(['Barcodes are being stretched between ' num2str(strMin) ' and ' num2str(strMax)]));

    commonLength = ceil(commonLength);

    import CBT.Consensus.Core.convert_barcodes_to_common_length;
    import CBT.Consensus.Core.convert_bitmasks_to_common_length;

    for i=1:numBar % change this to a simpler function
           % here interpolate both barcode and bitmask 
        lenBarTested = length(barcodeGenMol{i}.rawBarcode);
        barcodeGenMol{i}.stretchedBarcode = interp1(barcodeGenMol{i}.rawBarcode, linspace(1,lenBarTested,commonLength));
        barcodeGenMol{i}.rawBitmask = double(barcodeGenMol{i}.rawBitmask);
        barcodeGenMol{i}.rawBitmask(barcodeGenMol{i}.rawBitmask==0)=nan;
        barcodeGenMol{i}.stretchedrawBitmask = interp1(barcodeGenMol{i}.rawBitmask, linspace(1,lenBarTested,commonLength));
        barcodeGenMol{i}.stretchedrawBitmask(isnan(barcodeGenMol{i}.stretchedrawBitmask))=0;
        barcodeGenMol{i}.stretchedrawBitmask
%         [barcodeGenMol{i}.stretchedBarcode] = cell2mat(convert_barcodes_to_common_length({barcodeGenMol{i}.rawBarcode}, commonLength));
%         barcodeGenMol{i}.rawBitmask(barcodeGenMol{i}.rawBitmask==0)=nan;
%         [barcodeGenMol{i}.stretchedrawBitmask] = cell2mat(convert_bitmasks_to_common_length({barcodeGenMol{i}.rawBitmask}, commonLength));
%         barcodeGenMol{i}.stretchFactor = stretchings(i);
    end  
end
  
    
% now these we can output either as kymo, or ? 

% now, can run consensus for these barcodes
import CBT.Hca.Core.gen_consensus;
consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGenMol,sets);


import CBT.Hca.UI.Helper.select_all_consensus;
[consensusStructs.consensusStruct,sets] = select_all_consensus(consensusStructs,sets);

% consensusStructs.consensusStruct = consensusStruct;

% if we want to plot
import CBT.Hca.Export.plot_consensus_concentrically;
plot_consensus_concentrically(consensusStructs,barcodeGen)

%     
% import CBT.Hca.Core.gen_consensus_plot;
% gen_consensus_plot(consensusStructs,sets);

% take some index / or loop through all, might be slower
% idx = 8;

% 
% % we use kymoStructs to split each kymo to it's separate kymostruct, then
% % we can generate barcodes. We save each row of unaligned kymograph as
% % aligned kymo. Need to compute edges for the unaligned kymo (instead of
% % aligned)
% import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
% [ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(kymoStructs{idx}.unalignedKymo), sets.edgeDetectionSettings);
%      
% kymoSingleFrameStruct = [];
% sz = size(kymoStructs{idx}.unalignedKymo);
% kymos =mat2cell(kymoStructs{idx}.unalignedKymo,ones(1,sz(1)),sz(2));
% 
% for i=1:size(kymoStructs{idx}.unalignedKymo,1)
%     kymoSingleFrameStruct{i}.alignedKymo = double(kymos{i});
%     kymoSingleFrameStruct{i}.name = kymoStructs{idx}.name;
%     kymoSingleFrameStruct{i}.leftEdgeIdxs = leftEdgeIdxs(i);
%     kymoSingleFrameStruct{i}.rightEdgeIdxs = rightEdgeIdxs(i);
% end
% % 
% 
% % now we can generate barcodes
% import CBT.Hca.Core.gen_barcodes;
% barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoSingleFrameStruct, sets);
%  
% % barcode stretching factors. Note that edges are calculated only
% % approximately, but some tendencies in the stretchings can be seen here.
% % i.e stretching not completely random in between rows 
% figure,plot(cellfun(@(x) x.stretchFactor,barcodeGen))


