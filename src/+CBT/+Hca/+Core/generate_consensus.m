function [ consensusStruct ] = generate_consensus( barcodes, bitmasks, bgMeanApprox, sets )
    % generate_consensus is a function which finds a consensus barcode
    % based on a molecule structure, defined in molStruct, and set of
    % parameters sets, given as an input. The outpus is the consensus
    % structure, given in consensusStruct

    % input: barcodes, bitmasks, bgMeanApprox, sets
    % output: consensusStruct

    % code by Albertas Dvirnas, edited 09/11/18
    
    for k=1:length(bgMeanApprox)
        if isnan(bgMeanApprox{k})
            bgMeanApprox{k} = 0;
        end
    end

    % number of barcodes
    numBarcodes = length(barcodes);

    % lengths of barcodes
    lengths = cellfun(@length,barcodes);

    % Future possibility: allow barcodes of different lengths as an input
    % (for example for experiments with different conditions)
    % Therefore do not force stretching. That should be supported by HCM
    % project
    
    % Based on the choice of normalization, normalize the barcodes
    % before putting them into the matrix
    if strcmp(sets.consensus.barcodeNormalization, 'background')
        barcodeNormalizationFunction = @(bc, bg) (bc - bg);
       % barcodeNormalizationFunction = @(x) cellfun(bcInnerFunc, x, rawBgs, 'UniformOutput', 0);
    elseif strcmp(sets.consensus.barcodeNormalization, 'bgmean')
        barcodeNormalizationFunction = @(bc, bg) ((bc - bg) / nanmean(bc - bg));
       % barcodeNormalizationFunction = @(x) cellfun(bcInnerFunc, x, rawBgs, 'UniformOutput', 0);
    elseif strcmp(sets.consensus.barcodeNormalization, 'zscore')
        barcodeNormalizationFunction = @(bc,bg) zscore(bc-bg);
	elseif strcmp(sets.consensus.barcodeNormalization, 'original')
        barcodeNormalizationFunction = @(bc,bg) bc;
    end

    % define barcode matrix and bitmask. The length is that of the length
    % of longest barcode. 
    rawBar = zeros(numBarcodes,max(lengths));
    rawBit = zeros(numBarcodes,max(lengths));
      
    % put all the barcodes in the matrix
    for j=1:numBarcodes
        rawBar(j,1:lengths(j)) = barcodeNormalizationFunction(barcodes{j},bgMeanApprox{j});
        rawBit(j,1:lengths(j)) = bitmasks{j};
        rawBar(j,logical(~rawBit(j,:))) = 0;
    end

    % the things that are needed for the barcode comparisons are the
    % maximal coefficient matrix maxcoef, the orientation matrix or, and
    % the positional shift matris pos. 
    maxcoef = zeros(numBarcodes-1,numBarcodes-1);
    or = zeros(numBarcodes-1,numBarcodes-1);
    pos = zeros(numBarcodes-1,numBarcodes-1);
    
	%import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
	import SignalRegistration.masked_pcc_corr;
   % G = digraph();

    % run from the first barcode to the almost last one
    for barcodeIdxA = 1:numBarcodes-1
        barcodeA = rawBar(barcodeIdxA,:);
        bitmaskA = rawBit(barcodeIdxA,:);
        for barcodeIdxB = barcodeIdxA+1:numBarcodes
            barcodeB = rawBar(barcodeIdxB,:);
            bitmaskB = rawBit(barcodeIdxB,:);
            % compute correlatiton coefficients
            xcorrs = masked_pcc_corr(barcodeA,barcodeB,bitmaskA,bitmaskB);
            %xcorrs = get_no_crop_lin_circ_xcorrs(barcodeA,barcodeB,bitmaskA,bitmaskB);

            % compute maximum, position and orientation for the best
            % position
            [f,s] = max(xcorrs);          
            [ b, ix ] = sort( f(:), 'descend' );
            maxcoef(barcodeIdxA,barcodeIdxB-1) = b(1);
            or(barcodeIdxA,barcodeIdxB-1) = s(ix(1));
            pos(barcodeIdxA,barcodeIdxB-1) = ix(1);
        end     
    end

    % define new  barcode matrices that will change over time  and will
    % store the averaged barcodes
	rawBar2 = rawBar;
	rawBit2 = rawBit;

    % This matrix will store the information on which barcodes are averaged
    % in which row of the matrix which was defined a step before
     barToAverage = eye(numBarcodes,numBarcodes);
    % this keeps track on which barcodes are still not joined with
    % something
     barInd =1:numBarcodes;

     origPos = zeros(size(rawBar));
     origPos(:,1) = 1;
     origPos(:,2) = 2;

     %
     consensusStruct.treeStruct.maxCorCoef = zeros(1,numBarcodes-1);

 
    % Here we go through the barcodes until the global cluster (out of all
    % the barcodes) is computed, i.e. until we have the whole tree. 
    % Tree cutting is can be done here or left for post-processing.  
     for bb=1:numBarcodes-1
         
     
         % find which barcodes should be merged
         [consensusStruct.treeStruct.maxCorCoef(bb), I] = max(maxcoef(:));

         % 
        [I_row, I_col] = ind2sub(size(maxcoef),I);

        consensusStruct.treeStruct.clusteredBar{bb}{1} = find(barToAverage(I_row,:));
        consensusStruct.treeStruct.clusteredBar{bb}{2} = find(barToAverage(I_col+1,:));

        %consensusStruct.treeStruct.clusteredBar{bb} = strcat([mat2str(find(barToAverage(I_row,:))) ',' mat2str(find(barToAverage(I_col+1,:)))]);

        % todo: check which one is averaged over more barcodes, and 
        % here make a decision which one is to be removed
        if sum(barToAverage(I_row,:)) < sum(barToAverage(I_col+1,:))
            rowToRemove = I_row;
            rowToKeep = I_col+1;
        else
            rowToRemove = I_col+1;
            rowToKeep = I_row;
        end
        
        % indices of the averaged barcode we remove
        [a] = find(barToAverage(rowToRemove,:));

        % all the barcodes that will be placed here
        barToAverage(rowToKeep,:) = max(barToAverage(rowToKeep,:),barToAverage(rowToRemove,:));
         
        % We need to align all the barcodes from the one that is being
        % removed
        % TODO: instead of always having I_col > I_row, choose what to
        % align based on which averaged barcode has more barcodes in it. 
        % This also would give some accuracy and speed ! 

         % go through all barcodes that were in removed average
        for i=1:length(a)
            if rowToRemove > rowToKeep
                shiftInd = pos(rowToKeep,rowToRemove-1)-1;

                % first shift second barcode to the correct place
                rawBar(a(i),:) = circshift(rawBar(a(i),:),[0,-shiftInd]);
                rawBit(a(i),:) = circshift(rawBit(a(i),:),[0,-shiftInd]);
                origPos(a(i),:) = circshift(origPos(a(i),:),[0,-shiftInd]);
                if or(rowToKeep, rowToRemove-1) == 2
                    rawBar(a(i),:) = fliplr(rawBar(a(i),:));
                    rawBit(a(i),:) = fliplr(rawBit(a(i),:));
                    origPos(a(i),:) = fliplr(origPos(a(i),:));
                end
            end

            if rowToRemove < rowToKeep   
                shiftInd = pos(rowToRemove,rowToKeep-1)-1;

                if or(rowToRemove,rowToKeep-1) == 2
                    rawBar(a(i),:) = fliplr(rawBar(a(i),:));
                    rawBit(a(i),:) = fliplr(rawBit(a(i),:));
                    origPos(a(i),:) = fliplr(origPos(a(i),:));
                end
                
                rawBar(a(i),:) = circshift(rawBar(a(i),:),[0,shiftInd]);
                rawBit(a(i),:) = circshift(rawBit(a(i),:),[0,shiftInd]);
                origPos(a(i),:) = circshift(origPos(a(i),:),[0,shiftInd]);
            end
        end

        consensusStruct.treeStruct.barMatrix{bb} = logical(barToAverage(rowToKeep,:));
                
        positions = origPos(logical(barToAverage(rowToKeep,:)),:);
        [~,stPos]=find(positions==1);
        [~,stOr]=find(positions==2);
        barOr = (stPos > stOr) + 1;
        consensusStruct.treeStruct.barOrientation{bb} = [stPos barOr];
        % now we substitute I_row with this new barcode
        barcodesToAverage = rawBar(consensusStruct.treeStruct.barMatrix{bb},:);
        barcodesToAverage(~rawBit(consensusStruct.treeStruct.barMatrix{bb},:)) = nan;
        rawBar2(rowToKeep,:) = nanmean(barcodesToAverage);
        consensusStruct.treeStruct.barcodes{bb} = barcodesToAverage;
        
        rawBit2(rowToKeep,:) = max(rawBit(logical(barToAverage(rowToKeep,:)),:));
        rawBar2(rowToKeep,logical(~rawBit2(rowToKeep,:))) = 0;
      
        % todo: add an if statement so we don't recompute for rowToRemove
        
         % recompute pccs for barcode with index lower than rowToKeep
        for barcodeIdxA = 1:rowToKeep-1
            barcodeA = rawBar2(barcodeIdxA,:);
            bitmaskA = rawBit2(barcodeIdxA,:);
            barcodeB = rawBar2(rowToKeep,:);
            bitmaskB = rawBit2(rowToKeep,:);

            xcorrs = masked_pcc_corr(barcodeA,barcodeB,bitmaskA,bitmaskB);
            [f,s] = max(xcorrs);          
            [ b, ix ] = sort( f(:), 'descend' );
            indx = b(1) ;
            maxcoef(barcodeIdxA,rowToKeep-1) = indx;
            or(barcodeIdxA,rowToKeep-1) =s(ix(1));
            pos(barcodeIdxA,rowToKeep-1) = ix(1);
        end
        
         % recompute pccs for barcode with index higher than rowToKeep.
         % the total number of barcodes is numBarcodes-bb+1
        for barcodeIdxA = rowToKeep+1:numBarcodes-bb+1
            barcodeA = rawBar2(rowToKeep,:);
            bitmaskA = rawBit2(rowToKeep,:);
            barcodeB = rawBar2(barcodeIdxA,:);
            bitmaskB = rawBit2(barcodeIdxA,:);

            xcorrs = masked_pcc_corr(barcodeA,barcodeB,bitmaskA,bitmaskB);
            [f,s] = max(xcorrs);          
            [ b, ix ] = sort( f(:), 'descend' );
            indx = b(1) ;
            maxcoef(I_row,barcodeIdxA-1) = indx;
            or(I_row,barcodeIdxA-1) =s(ix(1));
            pos(I_row,barcodeIdxA-1) = ix(1);
        end
        
        % remove the extra row
        
        % remove rowToRemove
        rawBar2(rowToRemove,:) = [];
        rawBit2(rowToRemove,:) = [];
        
        % rowToRemove is now deleted.
        barToAverage(rowToRemove,:) = [];
        barInd(rowToRemove) = [];

        
        % In case I_col+1 is not the last column, we remove it from these..
        % this should always be satisfied?
        if rowToRemove <= size(maxcoef,1)
            maxcoef(rowToRemove,:) = [];
            pos(rowToRemove,:) = [];
            or(rowToRemove,:) = [];
        end

        % If I_col is not the first row, we remove as well
        % (Note it would change a bit if I_col or I_row could be removed -
        % now only I_col can be removed)
        try
            maxcoef(:,rowToRemove-1) = [];
            pos(:,rowToRemove-1) = [];
            or(:,rowToRemove-1) = [];
        catch % if first barcode is removed, we remove 2nd
            maxcoef(:,rowToRemove) = [];
            pos(:,rowToRemove) = [];
            or(:,rowToRemove) = [];
        end
        
     %    consensusStruct.treeStruct.bitMat{bb}= rawBit;
      %  consensusStruct.treeStruct.barMat{bb}= rawBar;
     end

    consensusStruct.time = datetime;
    
end

