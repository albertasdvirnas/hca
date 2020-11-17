function [] = pregen_zero_model_fft_from_prompted_barcodes()
    [barcodeFilenames, dirpath] = uigetfile('*.txt', ...
      'Select barcodes', 'MultiSelect', 'on');
    aborted = isequal(dirpath, 0);
    if aborted
        return;
    end
    if not(iscell(barcodeFilenames))
        barcodeFilenames = {barcodeFilenames};
    end
    barcodeFilepaths = fullfile(dirpath, barcodeFilenames(:));
    zeroModelBarcodes = cellfun(@importdata, barcodeFilepaths, 'un', 0);
    
    returnValue = inputdlg( ...
      'Maximum expected barcode length', ...
      'Select maximum expected barcode length');
    if isempty(returnValue)
      error('Unexpected returned value from prompt.');
    else
      maxExpectedBarcodeLength = str2double(returnValue{1});
    end
    zeroModelBarcodes = cellfun(@(x) x(cell2mat(arrayfun(@(i) ...
      i-maxExpectedBarcodeLength+1:i, ...
      maxExpectedBarcodeLength:maxExpectedBarcodeLength:length(x), ...
      'un', 0)')), ...
      zeroModelBarcodes, 'un', 0);
    zeroModelBarcodes = vertcat(zeroModelBarcodes{:});
    zeroModelBarcodes = mat2cell(zeroModelBarcodes, ...
      ones(1, size(zeroModelBarcodes, 1)));
      
    meanBpExt_pixels = nan;

    import CBT.UI.prompt_should_rescale;
    shouldRescale = prompt_should_rescale();
    if shouldRescale
        zeroModelBarcodes = cellfun(@zscore, zeroModelBarcodes, 'UniformOutput', false);
    end

    fprintf('Generating a zero-model fft...\n');
    numBarcodes = length(zeroModelBarcodes);
    fftMags = cellfun(@(x) abs(fft(x)), zeroModelBarcodes, 'un', 0);
    fftMagsMat = cell2mat(fftMags);
    meanZeroModelFftFreqMags = sqrt(sum(fftMagsMat.^2, 1)/numBarcodes);

    import CBT.RandBarcodeGen.PhaseRandomization.export_fft_file;
    export_fft_file(meanZeroModelFftFreqMags, meanBpExt_pixels);

    fprintf('Finished generating ZM fft from barcodes\n')
end