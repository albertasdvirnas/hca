function fftPhaseRandFactors = generate_fft_phase_randomization_factors(fftLen, numRandomizations)
    if nargin < 2
        numRandomizations = 1;
    end
    randLen = floor((fftLen - 1)/2);
    imagTau = (2 * pi) * (1i);
    tmpRand = rand(randLen, numRandomizations);
    fftPhaseRandFactors1 = exp(imagTau * tmpRand);
    fftPhaseRandFactors2 = exp(-imagTau * flipud(tmpRand));
    hasEvenFftLen = (mod(fftLen, 2) == 0);
    onesRowVect = ones(1, numRandomizations);
    if hasEvenFftLen
        fftPhaseRandFactors = [onesRowVect; fftPhaseRandFactors1; onesRowVect; fftPhaseRandFactors2];
    else
        fftPhaseRandFactors = [onesRowVect; fftPhaseRandFactors1; fftPhaseRandFactors2];
    end
end