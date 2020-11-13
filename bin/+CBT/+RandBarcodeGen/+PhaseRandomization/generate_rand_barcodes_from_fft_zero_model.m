function phaseRandomizedBarcodes = generate_rand_barcodes_from_fft_zero_model(meanZeroModelFftFreqMags, numRandomizations, barcodeLen)
    if nargin < 2
        numRandomizations = 1;
    end
    if nargin < 3
        barcodeLen = [];
    end
    
    if not(isempty(barcodeLen))
        import CBT.RandBarcodeGen.PhaseRandomization.resize_fft_freq_mags;
        meanZeroModelFftFreqMags = resize_fft_freq_mags(meanZeroModelFftFreqMags, barcodeLen);
    end
    
    import CBT.RandBarcodeGen.PhaseRandomization.generate_fft_phase_randomization_factors;
    
    outputBarcodesSz = size(meanZeroModelFftFreqMags);
    fftLen = length(meanZeroModelFftFreqMags);
    
    fftPhaseRandFactors = generate_fft_phase_randomization_factors(fftLen, numRandomizations); %fftLen X numRand matrix (each randomization gets a column)
    fftZeroModelPhaseRandomized = repmat(meanZeroModelFftFreqMags(:), [1, numRandomizations]).*fftPhaseRandFactors;
    
    phaseRandomizedBarcodes = cell(numRandomizations, 1);
    for randomizationNum = 1:numRandomizations
        phaseRandomizedBarcodes{randomizationNum} = reshape(ifft(fftZeroModelPhaseRandomized(:, randomizationNum)), outputBarcodesSz);
    end
end