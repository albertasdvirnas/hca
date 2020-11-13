function [] = export_fft_file(meanZeroModelFftFreqMags, meanBpExt_pixels, outputMatFilepath)
    import CBT.RandBarcodeGen.PhaseRandomization.try_prompt_pregen_fft_mat_output_filepath;
    if nargin < 3
        [~, outputMatFilepath] = try_prompt_pregen_fft_mat_output_filepath();
    end
    if isempty(outputMatFilepath)
        return;
    end

    % Saves the FFT and the kbpPerPixel value
    meanFFT = meanZeroModelFftFreqMags; %#ok<NASGU>
    kbpPerPixel = round(1/meanBpExt_pixels)/1000; %#ok<NASGU>
    save(outputMatFilepath, 'meanFFT', 'kbpPerPixel');
end