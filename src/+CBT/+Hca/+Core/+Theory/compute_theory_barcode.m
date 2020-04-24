function [ theory, header] = compute_theory_barcode( name,sets)
    % compute_hca_theory_barcode Computes theory barcode for hca
    %
    %   Args:
    %        chr1 : memory mapped theory sequence
    %        sets : settings, overlapLength
    %   Returns:
    %        theoryCurveUnscaled_pxRes, bitmask, probSeq, theorSeq
    
    % in case only theory is provided without settings, we list here all
    % the settings instead of hardcoding them later in the code.
    if nargin < 2
        sets.theoryGen.psfSigmaWidth_nm = 300;   % psf width in nm
        sets.theoryGen.meanBpExt_nm = 0.3;       % mean base-pair extension
        sets.theoryGen.k = 200000;                         %length of the small fragment to do fft on
        sets.theoryGen.m = 26000;                          % number of base-pairs to leave at the edges
        sets.theoryGen.isLinearTF;               % is theory circular
    end

    % create an easier accessible file
    import CBT.Hca.Core.Theory.create_memory_struct;
    [chr1,header] = create_memory_struct(name);
    
    if length(chr1.Data) < sets.theoryGen.k
        disp(strcat(['Theory ' name ' skipped because length < theoryGen.k =' num2str(sets.theoryGen.k) ]));
        theory = [];
%         header = [];
        delete(chr1.Filename);
        clear chr1;
        return;
    end
        

  
	% use reproducible random numbers. Fix them for the whole calculation,
	% so that we wouldn't have strange periodic structures in the barcode
    %

    % psf in basepairs, we need to convolve with such Gaussian
    psfSigmaWidth_bps = sets.theoryGen.psfSigmaWidth_nm / sets.theoryGen.meanBpExt_nm;

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
    ts(ts > 4) = randi(4,1,sum(ts > 4));
    rng(s);
    
    
    % length of ts, after adding the enge of sequence effects of length
    % m/2 at the begining and the end of the sequence
	n = length(ts);
    nSeq = length(chr1.Data); % theory + extra things on the left and right

    % clear chr1
    delete(chr1.Filename);
    clear chr1;
    
    % number of pixels. This comes from the size of original sequence
    numPx = round(nSeq/pxSize);
    
    % initialize theory barcode
    theory = zeros(1,numPx);
    
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
    
%     cutPoints = [1 2*floor(pxSize):floor(pxSize):nSeq]; % how best to deal with this. Don't want to have different
        
    % current index for px theory
    curIdx = 1;
        
    % initialize kernel (of length k)
    y = zeros(k,1);
    % kernel values
    y(1:m) = images.internal.createGaussianKernel(psfSigmaWidth_bps, m);

    %    Reverse query: not needed if kernel is symetric!
    %     y = ker(end:-1:1); %Reverse the query
    %     y(m+1:k) = 0; %append zeros

    % Fourier transform of the kernel. 
    Y = fft(y);
        
    % Only save the pixelated theory values in vector "theory"

    leftOverTheory = [];
    
    % wrapper for different possible theories
    import CBT.Hca.Core.Theory.compute_theory_wrapper;

    % main loop to compute the theories in batches of length k.
    for j = 1:k-m:n-k+1 % this is k-m+1 when computing the correlation. Maybe want to keep the same notation, and have k-m+1 here too?
        % compute theory for a fragment, which is 
        % take k elements and
        x = compute_theory_wrapper(ts(j:j+k-1), sets);
        
        %The main trick of getting dot products in O(n log n) time
        % compute fft of x
        X = fft(x);
        
        % product
        Z = X.*Y;
        
        % inverser fourier
        z = ifft(Z);
        
        % now z(m:k) contains the theory elements j:j+k-m
        % but we don't want to save these
%         theory(1,j:j+k-m) = z(m:k);    
        % so instead we save tempTheory
        tempTheory = [leftOverTheory; z(m:k-1)]; 
        
        % Need to investigate a bit more carefully which index to take here. If we add +1 here,
        % then we also always add+1 in the top sum, but this shift the
        % indexes by one bp over time. Does that affect the final output in
        % any meaningful way? We want temp theories to be of length k-m!
        % every step check how many pixels from thyCurve_pxRes we can
        % compute. we have j:j+k-m
        % j+k-m-1=length(tempTheory) ?
        
        pts = zeros(1,numPx);
        pts(curIdx:end) = j+k-m-1 >= cutPointsR(curIdx:end);
        vals = find(pts);

        % now all the theory values are available, later we will only
        % keep the ones that are needed for the calculation
        for idd = 1:length(vals)
%             data = tempTheory(cutPointsL(vals(idd))-j+1+length(leftOverTheory):cutPointsR(vals(idd))-j+1+length(leftOverTheory));
            theory(curIdx) = mean(tempTheory(cutPointsL(vals(idd))-j+1+length(leftOverTheory):cutPointsR(vals(idd))-j+1+length(leftOverTheory)));
            curIdx = curIdx +1;
        end
        
        % for a pixel 1, one needs 1:pxSize, pxSize+1:2pxSize,
        % 2pxSize+1:3pxSize
        if curIdx <= length(theory) 
            leftOverTheory = tempTheory(cutPointsL(curIdx)-j+1+length(leftOverTheory):end);
        end

    end
        if isempty(j)
            j = 0; % if nothing was computed
            k = n;
        else
            % at the end there are some points left
            j = j+k-m-1;
            k = n-j; % number of points left
        end

        if k >= m % if k < m, there are not enough points on theory to compute more px values
            
            x = compute_theory_wrapper(ts(j+1:n), sets);

            %The main trick of getting dot products in O(n log n) time
            X = fft(x);

            y(k+1:end)= [];

            Y = fft(y);

            Z = X.*Y;
            z = ifft(Z);
            
            tempTheory = [leftOverTheory; z(m:k-1)];

          % All the rest should be one
            pts = zeros(1,length(cutPointsR));
            pts(curIdx:end) = 1;
            vals = find(pts);

            % now all the theory values are available, later we will only
            % keep the ones that are needed for the calculation
            for idd = 1:length(vals) 
                theory(curIdx) = mean(tempTheory(cutPointsL(vals(idd))-j+length(leftOverTheory):cutPointsR(vals(idd))-j+length(leftOverTheory)));
                curIdx = curIdx +1;
            end
        
%             leftOverTheory = tempTheory(cutPointsR(vals(idd))-j+1:end);   
        end
        
end

