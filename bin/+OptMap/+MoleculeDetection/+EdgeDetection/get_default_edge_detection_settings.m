function [edgeDetectionSettings] = get_default_edge_detection_settings(skipDoubleTanhAdjustment)
    edgeDetectionSettings.skipDoubleTanhAdjustment = skipDoubleTanhAdjustment;

    % TODO: stop hardcoding these things
    
    otsuApproxSettings.globalThreshTF = true;
    otsuApproxSettings.smoothingWindowLen = 5;
    otsuApproxSettings.imcloseHalfGapLen = 10; % default 5
    otsuApproxSettings.numThresholds = 1;
    otsuApproxSettings.minNumThresholdsFgShouldPass = 1;
    edgeDetectionSettings.otsuApproxSettings = otsuApproxSettings;

    % TODO Note (from Saair):
    %  stop hardcoding these things:
    %  D & E in fit for
    %  f(x) = A + F *(tanh( (x - B) * D ) - tanh( (x - C) * E ))
    %  are being constrained between these values
    %  with initial estimate as given by
    %   tanhStretchFactorInitGuess = 1.0;
    %   minTanhStretchFactor = 0.5 * tanhStretchFactorInitGuess;
    %   maxTanhStretchFactor = 1.5 * tanhStretchFactorInitGuess;
    %  but these should probably be parameters that are calculated
    %  dynamically using setting information such as the
    %  point spread function's sigma width and maybe the expected
    %  fluctuations in the  molecule's movement over the duration
    %  of the exposure time for kymograph timeframe intensity
    %  profiles or some values based on numerical analysis
    %  Also Note that if the curve is rescaled away from the 
    %  current resolution (pixels at some constant resolution of
    %  pixel lengths in nm), these values ought to be adjusted 
    %  even if there is no change in wavelength/PSF width/velocity/
    %  stillness/exposure time, but currently that isn't accounted
    %  for
    tanhStretchFactorInitGuess = 1.0;
    minTanhStretchFactor = 0.5 * tanhStretchFactorInitGuess;
    maxTanhStretchFactor = 1.5 * tanhStretchFactorInitGuess;
    tanhSettings.tanhStretchFactorInitGuess = tanhStretchFactorInitGuess; 
    tanhSettings.minTanhStretchFactor = minTanhStretchFactor;
    tanhSettings.maxTanhStretchFactor = maxTanhStretchFactor;
    edgeDetectionSettings.tanhSettings = tanhSettings;
end