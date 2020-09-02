function [ kymoStructs,barcodeGenMol,consensusStructs ] = align_kymos_consensus( sets, kymoStructs )
    % align_kymos
    % Runs alignment of kymographs. Currently two choices of
    % methods - ssdalign and nralign. Add possibility of more methods in  
    %
    %     Args:
    %         sets: settings structure
    %         unalignedKymos: unaligned kymographs
    % 
    %     Returns:
    %         alignedKymo: aligned kymographs
    %         leftEdgeIdxs: left edge indices of the molecule
    %         rightEdgeIdxs: left edge indices of the molecule
   
    disp('Starting kymo alignment...')
    tic %

    % fix num of timeframes?
    
    barcodeGenMol = cell(1,length(kymoStructs));
    
    % first compute aligned kymographs
    for idx=1:length(kymoStructs)

        import CBT.Hca.Core.Helping.create_kymo_struct;
        [kymoSingleFrameStruct] = create_kymo_struct(kymoStructs{idx},sets);   % now we can generate barcodes

        % for kymographs, when averaging use original intensities
        sets.consensus.barcodeNormalization= 'original';

        import CBT.Hca.Core.gen_barcodes;
        barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoSingleFrameStruct, sets);

        import CBT.Hca.Core.gen_consensus;
        consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);


        import CBT.Hca.UI.Helper.select_all_consensus;
        [consensusStructs.consensusStruct ,sets] = select_all_consensus(consensusStructs,sets);
        
        kymoStructs{idx}.consensusStructs = consensusStructs;
        kymoStructs{idx}.barcodeGen = barcodeGen;
        kymoStructs{idx}.alignedKymo = consensusStructs.treeStruct.barcodes{end};
        kymoStructs{idx}.barcode = consensusStructs.consensusStruct{end}.rawBarcode;
        kymoStructs{idx}.bitmask = consensusStructs.consensusStruct{end}.rawBitmask;
        barcodeGenMol{idx}.bgMeanApprox = mean(cellfun(@(x) x.bgMeanApprox,barcodeGen));
        barcodeGenMol{idx}.bgStdApprox =  mean(cellfun(@(x) x.bgStdApprox,barcodeGen));
        barcodeGenMol{idx}.rawBarcode = consensusStructs.consensusStruct{end}.rawBarcode;
        barcodeGenMol{idx}.rawBitmask = consensusStructs.consensusStruct{end}.rawBitmask;
        barcodeGenMol{idx}.name = barcodeGen{1}.name;
%   
    end

    

    numBar = length(barcodeGenMol);
    sets.consensus.barcodeNormalization= 'bgmean';

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
        end  
    end
  

    % now these we can output either as kymo, or ? 

    % now, can run consensus for these barcodes
    import CBT.Hca.Core.gen_consensus;
    consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGenMol,sets);


    import CBT.Hca.UI.Helper.select_all_consensus;
    [consensusStructs.consensusStruct,sets] = select_all_consensus(consensusStructs,sets);



%     consensusStructs.consensusStruct = consensusStruct;
    
    timePassed = toc;
    disp(strcat(['All kymos were aligned in ' num2str(timePassed) ' seconds']));

	%assignin('base','hcaSessionStruct',hcaSessionStruct)

end

