function [meanFftFreqMags] = gen_mean_fft_freq_mags(barcodes)
    % Generates mean fft frequency magnitudes from the ZM barcodes
    % This is only used within the main program since no
    %  kbp/pixel calculation is done
    import CBT.RandBarcodeGen.PhaseRandomization.resize_fft_freq_mags;
    
    % Find longest ZM barcode
    barcodeLens = cellfun(@length, barcodes);
    newFftLen = max(barcodeLens);

    % Get the mean FT of the reference barcodes
    numBarcodes = length(barcodes);
    resizedBarcodeFftFreqMagsMat = zeros(numBarcodes, newFftLen);
    for barcodeNum = 1:numBarcodes
        currBarcode = barcodes{barcodeNum};
        oldFftFreqMags = abs(fft(currBarcode));
        resizedFftFreqMags = resize_fft_freq_mags(oldFftFreqMags, newFftLen);
        
        resizedBarcodeFftFreqMagsMat(barcodeNum, :) = resizedFftFreqMags;
    end
    meanFftFreqMags = sqrt(sum(resizedBarcodeFftFreqMagsMat .^ 2, 1) / numBarcodes);
end