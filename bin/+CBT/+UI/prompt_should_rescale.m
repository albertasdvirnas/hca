function [shouldRescale] = prompt_should_rescale()
    % Choose which kind of length input
    optShouldRescale = 'Yes';
    optShouldntRescale = 'No';
    optDefault = optShouldRescale;
    shouldRescalePromptMsg = 'Rescale the ZM barcodes? Press "No" only if the barcodes already are rescaled and you want to save time.';
    rescaleChoice = questdlg(shouldRescalePromptMsg, 'Rescale?', optShouldRescale, optShouldntRescale, optDefault);
    shouldRescale = strcmp(optShouldRescale, rescaleChoice);
end