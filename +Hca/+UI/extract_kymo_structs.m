function [ kymoStructs ] = extract_kymo_structs(rawKymos, rawKymoFilepaths, pixelsWidths_bps,alignedKymos,alignedKymosStretchFactors,shiftAlignedKymos)
    if nargin < 4
        passesFilters = ones(length(rawKymos),1);
        alignedKymos = cell(length(rawKymos),1);
        alignedKymosStretchFactors = cell(length(rawKymos),1);
        shiftAlignedKymos = cell(length(rawKymos),1);
    end
    
    dataType = 'Kymo';
    import Fancy.Utils.data_hash;
    dataHashes = cellfun(@data_hash, rawKymos, 'UniformOutput', false);
    %displayNames = strcat(srcFilenames, arrayfun(@(fileMoleculeIdx) [' Mol #', num2str(fileMoleculeIdx)], rawKymoFileMoleculeIdxs, 'UniformOutput', false));

    kymoStructs = struct(...
        'passesFilters', passesFilters,...
        'filePath', rawKymoFilepaths,...
        'unalignedKymo', rawKymos,...
        'dataHash', dataHashes,...
        'alignedKymos', alignedKymos,...
        'alignedKymosStretchFactors', alignedKymosStretchFactors,...
        'shiftAlignedKymos', shiftAlignedKymos,...
        'bpsPerPx', pixelsWidths_bps,...
        'type', dataType...
    );
    kymoStructs = arrayfun(@(kymoStruct) kymoStruct, kymoStructs, 'UniformOutput', false);
end
