function  [xcorrs] = unmasked_spearman_fft(A,B,W1)
%UNMASKED_SPEARMAN_FFT Summary of this function goes here
%   Detailed explanation goes here

% draft .. UNTESTED..

    Aw = A(logical(W1));
    % now convert Aw to rank
    AwRank = zscore(compute_ranks(Aw)); % or zscore(,1);
    BRank = compute_ranks(B); % or zscore(,1);

    Awrev = fliplr(AwRank);
    
%     shortVecCut = zscore(shortVec(logical(shortVecBit)));

    shortLength = size(AwRank,2);
    longLength = size(B,2);
    
    % here ranks of B are fixed instead of being in windows..
    longfft = fft(BRank);
    conVec = conj(fft(Aw,longLength));
    % Forward cross correlations
    ccForward = (ifft(longfft.*conVec))/(shortLength-1);

    conVec = conj(fft(Awrev,longLength));
    % Backward cross correlations
    ccBackward = (ifft(longfft.*conVec))/(shortLength-1);
    
    conVec = conj(fft(ones(1,shortLength),longLength));
    movMean =  (ifft(longfft.*conVec))./shortLength;
    movStd =  (ifft(fft(longVec.^2).*conVec));

    stdForward = sqrt((movStd-shortLength.*movMean.^2)./(shortLength-1 ));

    ccForward = ccForward ./ stdForward; 
    ccBackward = ccBackward ./stdForward; % std is the same for both forward and backward case

    xx2 = circshift(ccBackward,[0,1-find(shortVecBit,1,'first')]);
    xx1 = circshift(ccForward,[0,1-find(shortVecBit,1,'first')]);
    xcorrs = [xx1;xx2];
end

