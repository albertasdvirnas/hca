function [moleculeStartEdgeIdxsApprox, moleculeEndEdgeIdxsApprox, mainKymoMoleculeMaskApprox] = basic_otsu_approx_main_kymo_molecule_edges(kymo, globalThreshTF, smoothingWindowLen, imcloseHalfGapLen, numThresholds, minNumThresholdsFgShouldPass)
    % BASIC_OTSU_APPROX_MAIN_KYMO_MOLECULE_EDGES - Attempts to find the start and
    %  end indices for the main molecule in the kymograph
    %
    %  Uses Otsu's method,
    %  some morphological operations, and component analysis,
    %  to try to separate the foreground from the background and find the
    %  "main" (i.e. largest) contiguous foreground component in each row
    %  (i.e. time frame) represented in the provided kymograph
    %
    % Inputs:
    %   kymo
    %   globalThreshTF (optional, defaults to false)
    %   smoothingWindowLen (optional, defaults to 1)
    %   imcloseHalfGapLen (optional, defaults to 0)
    %   numThresholds (optional, defaults to 1)
    %   minNumThresholdsFgShouldPass (optional, defaults to 1)
    %   
    % Outputs:
    %   moleculeStartEdgeIdxsApprox
    %   moleculeEndEdgeIdxsApprox
    %   mainKymoMoleculeMaskApprox
    %
    % Authors:
    %   Saair Quaderi
    %     (refactoring)
    %   Charleston Noble
    %     (original, algorithm)

    if nargin < 2
        globalThreshTF = false;
    end
    
    if nargin < 3
        smoothingWindowLen = 1;
    end

    if nargin < 4
        imcloseHalfGapLen = 0;
    end

    if nargin < 5
        numThresholds = 1;
    end

    if nargin < 6
        minNumThresholdsFgShouldPass = 1;
    end


    kymoSz = size(kymo);
    numFrames = kymoSz(1);
    numCols = kymoSz(2);
    mainKymoMoleculeMaskApprox = false(kymoSz);
    kymoSmooth = kymo;
    if smoothingWindowLen > 1
        kymoSmooth = conv2(kymoSmooth, ones(1, smoothingWindowLen)/smoothingWindowLen, 'same');
    end
    if globalThreshTF
        thresholdsArr = multithresh(kymoSmooth(:), numThresholds);
        fgMask = kymoSmooth >= thresholdsArr(minNumThresholdsFgShouldPass);
    else
        fgMask = false(size(kymoSmooth));
        for frameNum = 1:numFrames
            kymoSmoothRow = kymoSmooth(frameNum, :);
            thresholdsArr = multithresh(kymoSmoothRow, numThresholds);
            fgMask(frameNum, :) = kymoSmooth(frameNum, :) >= thresholdsArr(minNumThresholdsFgShouldPass);
        end
    end
    
    ccFg = bwconncomp(fgMask);
    [rowIdxLists, colIdxLists] = cellfun(@(pixelIdxList) ind2sub(ccFg.ImageSize, pixelIdxList), ccFg.PixelIdxList, 'UniformOutput', false);
    
    touchesEdgeMask = cellfun(@(colIdxList) any(colIdxList == 1) | any(colIdxList == ccFg.ImageSize(2)), colIdxLists);
    if not(all(touchesEdgeMask)) % remove all connected components that touch edges (run out of field of view), unless everything touches edges
        ccFg.PixelIdxList = ccFg.PixelIdxList(~touchesEdgeMask);
        ccFg.NumObjects = sum(~touchesEdgeMask);
        fgMask = false(ccFg.ImageSize);
        for objectIdx = 1:ccFg.NumObjects
            fgMask(ccFg.PixelIdxList{objectIdx}) = true;
        end
    end
    
    % does not do anything?
    if imcloseHalfGapLen > 0
        imcloseNhood = true(1, 1 + 2*imcloseHalfGapLen);
        % Remove small gaps
        fgMask = imclose(fgMask, imcloseNhood);
        ccFg = bwconncomp(fgMask);
    end
    [~, tmp_so] = sort(cellfun(@length, ccFg.PixelIdxList), 'descend');
    ccFg.PixelIdxList = ccFg.PixelIdxList(tmp_so);
    [rowIdxLists, colIdxLists] = cellfun(@(pixelIdxList) ind2sub(ccFg.ImageSize, pixelIdxList), ccFg.PixelIdxList, 'UniformOutput', false);
    
    % Pick the largest molecule that includes the center column
    mainKymoMoleculeMaskApprox = false(ccFg.ImageSize);
    idx = 1;
    crossesCenterColIdxTF = false;
    while idx <= ccFg.NumObjects
        crossesCenterColIdxTF = (min(colIdxLists{idx}) <= numCols/2) && (max(colIdxLists{idx}) >= numCols/2);
        if crossesCenterColIdxTF
            mainKymoMoleculeMaskApprox(ccFg.PixelIdxList{idx}) = true;
            break;
        end
        idx = idx + 1;
    end
    if not(crossesCenterColIdxTF) && (ccFg.NumObjects > 0) % if nothing crosses center column, just pick the biggest molecule
        idx = 1;
        mainKymoMoleculeMaskApprox(ccFg.PixelIdxList{idx}) = true;
    end
    
    frameNums = (1:numFrames)';
    
    unemptyRows = arrayfun(@(frameNum) any(mainKymoMoleculeMaskApprox(frameNum, :)), frameNums);
    moleculeStartEdgeIdxsApprox = NaN(numFrames, 1);
    moleculeEndEdgeIdxsApprox = NaN(numFrames, 1);
    moleculeStartEdgeIdxsApprox(unemptyRows) = arrayfun(@(frameNum) find(mainKymoMoleculeMaskApprox(frameNum, :), 1, 'first'), frameNums(unemptyRows));
    moleculeEndEdgeIdxsApprox(unemptyRows) = arrayfun(@(frameNum) find(mainKymoMoleculeMaskApprox(frameNum, :), 1, 'last'), frameNums(unemptyRows));
end