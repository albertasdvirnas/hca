function [resizedFftFreqMags] = resize_fft_freq_mags(oldFftFreqMags, newFftLen)
    oldFftFreqMags = abs(oldFftFreqMags);
    oldFftLen = length(oldFftFreqMags);
    if (oldFftLen == newFftLen)
        resizedFftFreqMags = oldFftFreqMags;
        return;
    end
    newLenIsEven = (mod(newFftLen, 2) == 0);
    sharedInterpLen = floor((newFftLen - 1)/2);

    nearlyHalfOldFftLen = ceil(oldFftLen / 2) - 1; % largest integer strictly smaller than half the length
    slightlyMoreThanHalfOldFftLen = oldFftLen - nearlyHalfOldFftLen; % smallest integer strictly larger than half the length
    
    % Interpolate
    fftZeroFreqIdx = 1;
    interpPts = linspace(1, nearlyHalfOldFftLen, sharedInterpLen);
    fftMagsTempHalf = oldFftFreqMags(fftZeroFreqIdx + (1:nearlyHalfOldFftLen));
    fftMagsTempHalfStretched = interp1(fftMagsTempHalf, interpPts);

    % "Fold" negative frequences over positive
    if newLenIsEven
        resizedFftFreqMags = [...
            oldFftFreqMags(fftZeroFreqIdx), ...
            fftMagsTempHalfStretched, ...
            oldFftFreqMags(slightlyMoreThanHalfOldFftLen), ... % TODO: note from saair, this placement seems like it could be wrong (investigate)
            fliplr(fftMagsTempHalfStretched) ...
        ];
    else
        resizedFftFreqMags = [...
            oldFftFreqMags(fftZeroFreqIdx), ...
            fftMagsTempHalfStretched, ...
            fliplr(fftMagsTempHalfStretched)];
    end
    normFactor = sqrt((newFftLen * (newFftLen - 1)) / sum(resizedFftFreqMags .^ 2));
    resizedFftFreqMags = resizedFftFreqMags * normFactor;
end