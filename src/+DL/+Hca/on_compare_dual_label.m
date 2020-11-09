function [rezMaxM, bestBarStretch, bestLength] = on_compare_dual_label( ...
  barcodeGen, ...
  theoryStruct, ...
  stretchFactors, ...
  numPixelsAroundBestTheoryMask, ...
  sets)

% load theory barcode txt file.
formatSpec = '%f';
fileID = fopen(theoryStruct.filename,'r');
theorBar1 = transpose(fscanf(fileID,formatSpec));
fclose(fileID);
fileID = fopen(theoryStruct.filename2,'r');
theorBar2 = transpose(fscanf(fileID,formatSpec));
fclose(fileID);

% Temp. hardcoded parameters
% TODO: Load tabulated values into sets
theoryNumPixels = 3088286401/108.33*sets.theory.nmbp;
Rq = 1/0.0433;
NTrialsCb = 2*length(stretchFactors)/6.3;
MuDots = .46;
XiDots = .012;
EtaDots = -1.16e-6;
AlphaDots = .44;
BetaDots = 78;
NTrialsDots = 2*length(stretchFactors);

rezMaxM = cell(1, size(barcodeGen, 2));
bestBarStretch = zeros(1, size(barcodeGen, 2));
bestLength = zeros(1, size(barcodeGen, 2));
% for all the barcodes run
parfor idx=1:size(barcodeGen, 2)
  
  % scoreMax stores the  maximum coefficients
  scoreMax = zeros(3,length(stretchFactors));
  
  % rezMaz stores the results for one barcode
  rezMax = cell(1,length(stretchFactors));
  
  % length of this barcode
  lenBarTested = length(barcodeGen{1, idx}.rawBarcode);
  
  % run the loop for the stretch factors
  for j=1:length(stretchFactors)
    % here interpolate both barcode and bitmask
    v = linspace(1,lenBarTested,lenBarTested*stretchFactors(j));
    barDense = interp1(barcodeGen{1, idx}.rawBarcode, v);
    barSparse = interp1(barcodeGen{2, idx}.rawBarcode, v);
    
    import DL.Hca.qGauss_ev_cdf
    xDense = MASS_PCC( ...
      theorBar1, ...
      barDense, ...
      2^(4+nextpow2(length(barDense))));
    zDense = -norminv(1-qGauss_ev_cdf(xDense, ...
      Rq/lenBarTested, ...
      NTrialsCb*(theoryNumPixels-lenBarTested)));
    doIncreasedPrecision = isinf(zDense) & zDense > 0;
    zDense(doIncreasedPrecision) = -norminv(double(vpa(1)-qGauss_ev_cdf(xDense(doIncreasedPrecision), ...
      Rq/lenBarTested, ...
      NTrialsCb*(theoryNumPixels-lenBarTested), true)));
    
    import DL.Hca.trunc_normal_ev_cdf
    xSparse = MASS_DOT_CC( ...
      theorBar2, ...
      barSparse, ...
      2^(4+nextpow2(length(barSparse))));
    zSparse = -norminv(1-trunc_normal_ev_cdf(xSparse, ...
      norminv(linspace(0.001, 0.999, 100), MuDots, ...
      XiDots+EtaDots*lenBarTested), ...
      NTrialsDots*(theoryNumPixels-lenBarTested), ...
      1/sqrt(AlphaDots*lenBarTested+BetaDots), 0, 1));
    doIncreasedPrecision = isinf(zSparse) & zSparse > 0;
    zSparse(doIncreasedPrecision) = -norminv(double(vpa(1)-trunc_normal_ev_cdf(xSparse(doIncreasedPrecision), ...
      norminv(linspace(0.001, 0.999, 100), MuDots, ...
      XiDots+EtaDots*lenBarTested), ...
      NTrialsDots*(theoryNumPixels-lenBarTested), ...
      1/sqrt(AlphaDots*lenBarTested+BetaDots), 0, 1, true)));
    
    import CBT.Hca.UI.Helper.get_best_parameters;
    [rezMax{j}.maxcoef, ...
      rezMax{j}.pos, ...
      rezMax{j}.or] = get_best_parameters((zDense + zSparse)/sqrt(2), ...
      3, length(barDense), 1, numPixelsAroundBestTheoryMask);
    [rezMax{j}.maxcoefDense, ...
      rezMax{j}.posDense, ...
      rezMax{j}.orDense] = get_best_parameters(zDense, ...
      3, length(barDense), 1, numPixelsAroundBestTheoryMask);
    [rezMax{j}.maxcoefSparse, ...
      rezMax{j}.posSparse, ...
      rezMax{j}.orSparse] = get_best_parameters(zSparse, ...
      3, length(barSparse), 1, numPixelsAroundBestTheoryMask);
    
    rezMax{j}.indcoef = [ ...
      arrayfun(@(i) zDense(rezMax{j}.or(i), rezMax{j}.pos(i)), 1:3); ...
      arrayfun(@(i) zSparse(rezMax{j}.or(i), rezMax{j}.pos(i)), 1:3)];
    
    rezMax{j}.secondPos = 0;
    rezMax{j}.lengthMatch = 0;
    rezMax{j}.dist = 0;
    
    scoreMax(1, j) = rezMax{j}.maxcoef(1);
    scoreMax(2, j) = rezMax{j}.maxcoefDense(1);
    scoreMax(3, j) = rezMax{j}.maxcoefSparse(1);
  end
  
  % TODO: Save best stretch/length for each score type
  
  [~,b] = nanmax(scoreMax(1,:));
  rezMaxM{idx} = rezMax{b};
  [~,c] = nanmax(scoreMax(2,:));
  rezMaxM{idx}.maxcoefDense = rezMax{c}.maxcoefDense;
  rezMaxM{idx}.posDense = rezMax{c}.posDense;
  rezMaxM{idx}.orDense = rezMax{c}.orDense;
  [~,d] = nanmax(scoreMax(3,:));
  rezMaxM{idx}.maxcoefSparse = rezMax{d}.maxcoefSparse;
  rezMaxM{idx}.posSparse = rezMax{d}.posSparse;
  rezMaxM{idx}.orSparse = rezMax{d}.orSparse;
  bestBarStretch(idx) = stretchFactors(b);
  bestLength(idx) = round(lenBarTested*stretchFactors(b));
  %   end
  
%   disp(num2str([idx rezMaxM{idx}.maxcoef(1) rezMaxM{idx}.maxcoefDense(1) rezMaxM{idx}.maxcoefSparse(1)]))
end
