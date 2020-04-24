function [ consensus,sets] = select_all_consensus( consensusStructs, sets )
    % select_all_consensus - generates all possible consensus pairs
    
    % 	Args:
    %      consensusStructs, sets 
    %   Returns:
    %       consensus,sets
    
    % output consensusStructs
    
    
    if (sets.genConsensus==0) || isempty(consensusStructs)
        consensus = [];
        return;
    end
    consensus = cell(1,length(consensusStructs.treeStruct.maxCorCoef));
    
    % loop through consensus tree structure
    for idx=1:length(consensusStructs.treeStruct.maxCorCoef)
        barsAveraged = sort(cell2mat(consensusStructs.treeStruct.clusteredBar{idx}));
        
        consensus{idx}.rawBarcode = nanmean(consensusStructs.treeStruct.barcodes{idx});
        consensus{idx}.indices = barsAveraged;

        bitm = isnan(consensusStructs.treeStruct.barcodes{idx});
        consensus{idx}.indexWeights = sum(~bitm);
        consensus{idx}.rawBitmask = sum(~bitm)> 3*std(sum(bitm));
        consensus{idx}.rawBarcode(~consensus{idx}.rawBitmask) = nan;
        consensus{idx}.stdErrOfTheMean = nanstd(consensusStructs.treeStruct.barcodes{idx})./sqrt(consensus{idx}.indexWeights);
        consensus{idx}.time = consensusStructs.time;  
        consensus{idx}.idx = idx;
    end

end

