function [aborted, outputMatFilepath] = try_prompt_pregen_fft_mat_output_filepath()
    % Choose output file name
    [outputMatFilename, outputDirpath] = uiputfile('*.mat', 'Save Pregen FFT As');
    aborted = isequal(outputDirpath, 0);
    if aborted
        outputMatFilepath = '';
        return;
    end
    outputMatFilepath = fullfile(outputDirpath, outputMatFilename);
end