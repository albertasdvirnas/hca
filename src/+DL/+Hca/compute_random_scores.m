function [bestScore, scoreMean] = compute_random_scores( ...
  randLong, ...
  randShort, ...
  stretchMax, ...
  barcodeType)

stretchFactor = 1-stretchMax:0.01:1+stretchMax;
slidingScoreVec = nan(1, length(stretchFactor));
slidingScoreMeanVec = nan(1, length(stretchFactor));

for i=1:length(stretchFactor)
  v = linspace(1, length(randShort), round(length(randShort)*stretchFactor(i)));
  interpBarcode = interp1(randShort, v);
  switch barcodeType
    case 'cb'
      scoreMat = MASS_PCC(randLong, interpBarcode, 2^(4+nextpow2(length(interpBarcode))));
    case 'dots'
      scoreMat = MASS_DOT_CC(randLong, interpBarcode, 2^(4+nextpow2(length(interpBarcode))));
  end
  slidingScoreVec(i) = nanmax(scoreMat(:));
  slidingScoreMeanVec(i) = nanmean(scoreMat(:));
end

bestScore = nanmax(slidingScoreVec(:));
scoreMean = nanmean(slidingScoreMeanVec(:));