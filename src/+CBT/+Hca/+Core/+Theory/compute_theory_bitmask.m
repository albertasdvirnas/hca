function [ bitmask] = compute_theory_bitmask( chr1,sets)
    % compute_theory_bitmask Computes hca theory bitmask. 
    %   0 in a region that has nan's
    %
    %   Args:
    %        chr1 : memory mapped theory sequence
    %        sets : settings, overlapLength
    %   Returns:
    %        theoryCurveUnscaled_pxRes, bitmask, probSeq, theorSeq
    

    % This function converts bpRes to pxRes
    pxSize = sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm;

    % maybe no need to re-write these?    
    k = sets.theoryGen.k; % length of the small fragment to do fft on
    m = sets.theoryGen.m; % number of base-pairs to leave at the edges
    circular = ~sets.theoryGen.isLinearTF;
    
    % load series to a new variable. This could be later removed
    % since we can directly access values of the ts without putting it
    % to memory
    
    % keep same as in the old convert_bpres_to_pxres for reproducibility
    % in that function, we add round(s) and round(2*s) to the beginning and
    % and of the sequence after computation to make it circular.
    % Here we
    % deal with circularity by adding m/2 values to the beginning and the
    % end. (extraL+ round(s) and extraR+round(2*s)) 
    
    % this is how things are averaged ..
       
    extraL = m/2+round(pxSize)+1; % this +1 comes from convert_bpres_to_pxres, but there's no real argument why to use it
    extraR = m/2+round(2*pxSize); 
    ts = [ones(extraL,1); chr1.Data; ones(extraR,1)];
    if circular % so if circular, we add the theory data instead of simple ones
        ts(1:extraL) = chr1.Data(end-extraL+1:end);
        ts(end-extraR+1:end) =  chr1.Data(1:extraR);
    end
        
    % This is moved to outer structure avoid repetitive (periodic) things in the theory barcode.
    s = rng;
    rng(0,'twister');
    ts(ts > 4) = 0;
    rng(s);
    
    % length of ts, after adding the enge of sequence effects of length
    % m/2 at the begining and the end of the sequence
    nSeq = length(chr1.Data); % theory + extra things on the left and right

    % number of pixels. This comes from the size of original sequence
    numPx = round(nSeq/pxSize);
    
    % initialize theory barcode
    bitmask = zeros(1,numPx);
    
    % follow notation from convert_bpRes_to_pxRes
%     theory(floor(nSeq/pxSize)) = 0;
    % cutPoints. after reaching these points, the value in a pixel of
    % theory is computed. the first one starts from 1. We don't end exactly
    % and nSeq, so there's some information loss (almost up to a pixel) 
    % especially if the barcode is circular. 
 
    % find the cutpoints in theory at which should we compute the theory
    cutPointsR = ones(1,numPx);
    cutPointsL = ones(1,numPx);
    for i = 1:length(cutPointsR)
        cutPointsR(i) = floor((i*pxSize)+pxSize); % so we take [-pxSize:pxSize] around the value of interest
        cutPointsL(i) = floor((i*pxSize)-pxSize+1);
    end
    % note that  cutPointsL(i):cutPointsR(i) will not always be the same
    % and changes a little bit, because we take the floor value ( so can
    % differ in +-1 basepair)
    
    for j = 1:length(cutPointsL)
        bitmask(j) = sum(ts(cutPointsL(j):cutPointsR(j))==0);
    end
    
        
end

