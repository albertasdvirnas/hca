function [ barcodeGen ] = gen_barcodes( kymoStructs, sets )
    % gen_barcodes
    % generates barcodes from aligned kymograph data
    %
    %     Args:
    %         sets: settings structure
    %         kymoStructs: structure file for kymographs
    % 
    %     Returns:
    %         barcodeGen: barcode data
    
    % number of barcodes
    numBar = length(kymoStructs);
    
    % predefine cell for barcode structure
    barcodeGen = cell(1,numBar);
    
    disp('Starting generating barcodes...')
    tic
    % this computes barcode structure for nralign or other alignment methods 
    import CBT.Hca.Core.gen_barcode_data;
    
    % in case we need to filter kymos
   	import CBT.Hca.Core.filter_kymos;

    for i=1:numBar 
        % in case the kymo needs to be filtered, run this. This is not
        % saved for later, maybe we want to 
        kymoStructs{i}.alignedKymo = filter_kymos(kymoStructs{i}.alignedKymo,sets.filterSettings);
        % generate barcode data
        [barcodeGen{i}] = gen_barcode_data(kymoStructs{i}.alignedKymo,kymoStructs{i}.leftEdgeIdxs, kymoStructs{i}.rightEdgeIdxs,sets.skipEdgeDetection);
        barcodeGen{i}.name =  kymoStructs{i}.name;
    end

    % now compute information score
    %informationScores = cell(1,length(hcaSessionStruct.unalignedKymos));
    for i=1:numBar
        tempMean = mean(kymoStructs{i}.alignedKymo(:,barcodeGen{i}.lE:barcodeGen{i}.rE),2); % mean across x
        barcodeGen{i}.infoscore.mean = mean(tempMean); % total mean
        if size(tempMean,1) ==1
            barcodeGen{i}.infoscore.std =  std(kymoStructs{i}.alignedKymo(:,barcodeGen{i}.lE:barcodeGen{i}.rE)); % total std
        else
            barcodeGen{i}.infoscore.std =  std(tempMean); % total std
        end

        barcodeGen{i}.infoscore.score = barcodeGen{i}.infoscore.mean+3* barcodeGen{i}.infoscore.std;
    end
    
    %% Now define bitmasks
    % untrusted number of pixels, also will change if stretchfactor is
    % introduced
    untrPx = sets.bitmasking.untrustedPx;

	import CBT.Bitmasking.generate_zero_edged_bitmask_row;
    % add bitmasks. TODO: Also consider if there could be NAN's in the middle of
    % the molecule (due to, i.e., fragmentation)
    for i=1:numBar
        % add standard bitmask where pixels in the beginning and the end
        % are ignored
        barcodeGen{i}.rawBitmask = generate_zero_edged_bitmask_row(length(barcodeGen{i}.rawBarcode),round(untrPx));
        % in case some untrusted pixels in the middle (i.e. due to two
        % molecules being close to each other), bitmask this region
        barcodeGen{i}.rawBarcode(isnan(barcodeGen{i}.rawBarcode)) = false;
    end
    
	%% Prestretch barcodes to the same length
    % convert to common length, if chosen
    if  sets.genConsensus == 1
        allLengths = cellfun(@(x) length(x.rawBarcode),barcodeGen);
        commonLength = ceil(mean(allLengths));
        stretchings = commonLength./allLengths;
        strMin = min(stretchings);
        strMax =  max(stretchings);
        disp(strcat(['Barcodes are being stretched between ' num2str(strMin) ' and ' num2str(strMax)]));

        commonLength = ceil(commonLength);

        
        import CBT.Consensus.Core.convert_barcodes_to_common_length;
        import CBT.Consensus.Core.convert_bitmasks_to_common_length;

        for i=1:numBar % change this to a simpler function
            [barcodeGen{i}.stretchedBarcode] = cell2mat(convert_barcodes_to_common_length({barcodeGen{i}.rawBarcode}, commonLength));
            [barcodeGen{i}.stretchedrawBitmask] = cell2mat(convert_bitmasks_to_common_length({barcodeGen{i}.rawBitmask}, commonLength));
            barcodeGen{i}.stretchFactor = stretchings(i);
        end  
    end
  
    timePassed = toc;
    disp(strcat(['All barcodes generated in ' num2str(timePassed) ' seconds']));

end

