function [kymoSingleFrameStruct] = create_kymo_struct(kymoStructs,sets)
    % create_kymo_struct creates kymo struct for individual rows for
    % multiframe kymos
    
    %   Args: 
    %       kymoStructs : kymo structure
    %       sets : settings
    %
    %   Returns:
    %       kymoSingleFrameStruct : output structure

% take some index / or loop through all, might be slower
% for idx = 1:length(kymoStructs)

    % we use kymoStructs to split each kymo to it's separate kymostruct, then
    % we can generate barcodes. We save each row of unaligned kymograph as
    % aligned kymo. Need to compute edges for the unaligned kymo (instead of
    % aligned)
    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [ leftEdgeIdxs,rightEdgeIdxs,~] = approx_main_kymo_molecule_edges(double(kymoStructs.unalignedKymo), sets.edgeDetectionSettings);

    sz = size(kymoStructs.unalignedKymo);
    kymos =mat2cell(kymoStructs.unalignedKymo,ones(1,sz(1)),sz(2));

    if sets.timeFramesNr == 0
        minV = inf;
    else
        minV = sets.timeFramesNr;
    end
    
    kymoSingleFrameStruct = cell(1, min(minV,size(kymoStructs.unalignedKymo,1)));

    for i=1:min(minV,size(kymoStructs.unalignedKymo,1))
        kymoSingleFrameStruct{i}.alignedKymo = double(kymos{i});
        kymoSingleFrameStruct{i}.name = kymoStructs.name;
        kymoSingleFrameStruct{i}.leftEdgeIdxs = leftEdgeIdxs(i);
        kymoSingleFrameStruct{i}.rightEdgeIdxs = rightEdgeIdxs(i);
    end
    % 
% end

