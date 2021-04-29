function rezMaxM = on_compare_dual_label(barcodeGen, theoryStruct, sets, externalAlignmentStruct)
  import CBT.Hca.Import.load_pval_struct;
  import DL.Hca.beta_ev_cdf
  import SignalRegistration.masked_pcc_corr
  import DL.Hca.trunc_normal_ev_cdf
  import CBT.Hca.UI.Helper.get_best_parameters;
  
  if isfield(externalAlignmentStruct, 'molIds')
    doExternal = true;
  else
    doExternal = false;
  end

  %% Change to balance precision/speed
  digits(64) % Number of digits precision.
  xThreshDense = 5;
  xThreshSparse = 4;
  quenchPointsLow = linspace(0.001, 0.999, 31);
  quenchPointsHi = linspace(0.001, 0.999, 51);
  stretchTolerance = 2;
  positionTolerance = 5;

  %% load theory barcode txt file.
  formatSpec = '%f';
  fileID = fopen(theoryStruct.filename, 'r');
  theorBar1 = transpose(fscanf(fileID, formatSpec));
  fclose(fileID);
  fileID = fopen(theoryStruct.filename2, 'r');
  theorBar2 = transpose(fscanf(fileID, formatSpec));
  fclose(fileID);

  try
    fileID = fopen(strrep(theoryStruct.filename, 'barcode', 'bitmask'), 'r');
    theorUnkCount = transpose(fscanf(fileID, formatSpec));
    fclose(fileID);
  catch
    theorUnkCount = zeros(size(theorBar1));
  end
  theorBitmask = true(size(theorBar1));

  %% Load zero-model params
  theoryNumPixels = 3088286401 / sets.theory.pixelWidth_nm * sets.theory.nmbp;
  numStretchFactors = length(sets.theory.stretchFactors);
  psfPxRounded = round(sets.theory.psfSigmaWidth_nm / sets.theory.pixelWidth_nm, 4);
  cbParamsPath = fullfile(sets.duallabel.paramFolder, strcat(sets.duallabel.denseBarcodeType, "_zero_model_params.txt"));
  dotsParamsPath = fullfile(sets.duallabel.paramFolder, compose("dots_%s_zero_model_params.txt", sets.duallabel.sparsePattern));
  [cbPsfInd, cbParams] = load_pval_struct(cbParamsPath);
  [dotsPsfInd, dotParams] = load_pval_struct(dotsParamsPath);
  cbParamsInList = ismember(psfPxRounded, cbPsfInd);
  dotParamsInList = ismember(psfPxRounded, dotsPsfInd);

  try
    cbParams = cbParams{cbParamsInList};
    cbParamLambda = cbParams(1);
    cbParamExa = cbParams(2);
    cbParamNu = cbParams(3);
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", psfPxRounded, cbParamsPath));
  end

  try
    dotParams = dotParams{dotParamsInList};
    dotParamLambda = dotParams(1);
    dotParamExa = dotParams(2);
    dotParamMu = dotParams(3);
    dotParamXi = dotParams(4);
    dotParamSigma = dotParams(5);
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", psfPxRounded, dotsParamsPath));
  end

  %% Calculate scores
  numBarcodes = length(barcodeGen);
  rezMaxM = cell(1, numBarcodes);
  stretchFactors = sets.theory.stretchFactors;
  doBitmasking = sets.duallabel.doBitmask;
  numPxMatchDiff = sets.duallabel.numPxMatchDiff;
  barcodeGenDense = barcodeGen(1, :);
  barcodeGenSparse = barcodeGen(2, :);

  for idx = 1:numBarcodes

    rezMax = struct();

    xDenseAll = cell(1, numStretchFactors);
    xSparseAll = cell(1, numStretchFactors);
    zDenseAll = cell(1, numStretchFactors);
    zSparseAll = cell(1, numStretchFactors);

    denseBar = barcodeGenDense{idx}.rawBarcode;
    sparseBar = barcodeGenSparse{idx}.rawBarcode;
    bitmask = barcodeGenDense{idx}.rawBitmask;
    lenBarDense = length(denseBar);
    lenBarSparse = length(sparseBar);

    % run the loop for the stretch factors
    for j = 1:numStretchFactors
      
      % here interpolate both barcode and bitmask
      lenBarDenseStretched = round(lenBarDense * stretchFactors(j));
      lenBarSparseStretched = round(lenBarSparse * stretchFactors(j));
      vDense = linspace(1, lenBarDense, lenBarDenseStretched);
      vSparse = linspace(1, lenBarSparse, lenBarSparseStretched);
      thisBarDense = interp1(denseBar, vDense);
      thisBarSparse = interp1(sparseBar, vSparse);
      thisBitmask = bitmask(round(vDense));

      % Dense scores
      if doBitmasking
        xDense = masked_pcc_corr(thisBarDense, theorBar1, thisBitmask, theorBitmask);
        nuLenBar = sum(thisBitmask);
      else
        xDense = MASS_PCC(theorBar1, thisBarDense, 2^(4 + nextpow2(lenBarDenseStretched)));
        nuLenBar = lenBarDenseStretched;
      end

      lambdaEffDense = cbParamLambda * cbParamExa * numStretchFactors * (theoryNumPixels - lenBarDenseStretched);
      doZScore = xDense > nanstd(xDense(:)) * xThreshDense;
      zDense = zeros(size(xDense));
      zDense(doZScore) = -2 * log(vpa(1) - beta_ev_cdf(xDense(doZScore), max(4, cbParamNu * nuLenBar), 1, lambdaEffDense));
      
      doIncreasedPrecision = isinf(zDense) & zDense > 0;
      zDense(doIncreasedPrecision) = -2 * log(vpa(1) - beta_ev_cdf(xDense(doIncreasedPrecision), max(4, cbParamNu * nuLenBar), 1, lambdaEffDense, true));

      % Sparse scores
      xSparse = MASS_DOT_CC(theorBar2, thisBarSparse, 2^(4 + nextpow2(lenBarSparseStretched)));
      
      lambdaEffSparse = dotParamLambda * dotParamExa * numStretchFactors * (theoryNumPixels - lenBarSparseStretched);
      doZScore = xSparse > nanmean(xSparse(:)) + nanstd(xSparse(:)) * xThreshSparse;
      zSparse = zeros(size(xSparse));
      zSparse(doZScore) = -2 * log(vpa(1) - trunc_normal_ev_cdf(xSparse(doZScore), norminv(quenchPointsLow, dotParamMu, 1 / sqrt(dotParamXi * lenBarSparseStretched)), 1 / sqrt(dotParamSigma * lenBarSparseStretched), lambdaEffSparse, 0, 1));
      
      doIncreasedPrecision = isinf(zSparse) & zSparse > 0;
      zSparse(doIncreasedPrecision) = -2 * log(vpa(1) - trunc_normal_ev_cdf(xSparse(doIncreasedPrecision), norminv(quenchPointsHi, dotParamMu, 1 / sqrt(dotParamXi * lenBarSparseStretched)), 1 / sqrt(dotParamSigma * lenBarSparseStretched), lambdaEffSparse, 0, 1, true));

      xDenseAll{j} = xDense;
      xSparseAll{j} = xSparse;
      zDenseAll{j} = zDense;
      zSparseAll{j} = zSparse;

    end

    % Optimise combined score
    zDualAll = cell(1, numStretchFactors);
    posOffset = cell(1, numStretchFactors);
    stretchOffset = cell(1, numStretchFactors);

    for j = 1:numStretchFactors
      stretchIds = max(1, j - stretchTolerance):min(numStretchFactors, j + stretchTolerance);

      thisSparseZ = zSparseAll{j};
      thisDualZ = zeros([size(zSparseAll{j}), length(stretchIds), 2 * positionTolerance + 1]);

      for k = 1:length(stretchIds)
        for m = 1:2 * positionTolerance + 1
          shift = -(m - positionTolerance - 1);
          thisDualZ(:, :, k, m) = nansum(cat(3, thisSparseZ, circshift(zDenseAll{stretchIds(k)}, [0 shift])), 3);
        end
      end

      tmpScoreMat1 = nanmax(thisDualZ, [], 4);
      tmpScoreMat2 = nanmax(thisDualZ, [], 3);
      [~, tmpStretch] = nanmax(tmpScoreMat1, [], 3);
      [zDualAll{j}, tmpPos] = nanmax(tmpScoreMat2, [], 4);
      stretchOffset{j} = reshape(stretchIds(tmpStretch) - j, size(zSparseAll{j}));
      posOffset{j} = reshape(tmpPos - positionTolerance - 1, size(zSparseAll{j}));
      
      isUnchanged = zDualAll{j} == zDenseAll{j} | zDualAll{j} == zSparseAll{j};
      stretchOffset{j}(isUnchanged) = 0;
      posOffset{j}(isUnchanged) = 0;
    end

    % Get best stretchF for individual barcodes
    tmpLongestBar = round(max(lenBarDense, lenBarSparse) * stretchFactors(end));
    [~, b] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zDualAll));
    [~, c] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zDenseAll));
    [~, d] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zSparseAll));

    rezMax.dual = struct();
    rezMax.dual.bestBarStretch = stretchFactors(b);
    rezMax.dual.bestLength = round(max(lenBarDense, lenBarSparse) * stretchFactors(b));
    [rezMax.dual.maxcoef, rezMax.dual.pos, rezMax.dual.or] = get_best_parameters(zDualAll{b}, 3, rezMax.dual.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.dual.posOffset = arrayfun(@(i) posOffset{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3);
    rezMax.dual.stretchOffset = arrayfun(@(i) stretchOffset{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3);
    optPos = max(1, rezMax.dual.pos + rezMax.dual.posOffset);
    optStr = b + rezMax.dual.stretchOffset;
    rezMax.dual.optStr = stretchFactors(optStr);
    rezMax.dual.maxcoefParts = [arrayfun(@(i) zDenseAll{optStr(i)}(rezMax.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) zSparseAll{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3)];
    rezMax.dual.maxcoefPartsCC = [arrayfun(@(i) xDenseAll{optStr(i)}(rezMax.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) xSparseAll{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3)];

    rezMax.dense = struct();
    rezMax.dense.bestBarStretch = stretchFactors(c);
    rezMax.dense.bestLength = round(lenBarDense * stretchFactors(c));
    [rezMax.dense.maxcoef, rezMax.dense.pos, rezMax.dense.or] = get_best_parameters(zDenseAll{c}, 3, rezMax.dense.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.dense.maxcoefParts = [arrayfun(@(i) zDenseAll{c}(rezMax.dense.or(i), rezMax.dense.pos(i)), 1:3); arrayfun(@(i) zSparseAll{c}(rezMax.dense.or(i), rezMax.dense.pos(i)), 1:3)];
    rezMax.dense.maxcoefPartsCC = [arrayfun(@(i) xDenseAll{c}(rezMax.dense.or(i), rezMax.dense.pos(i)), 1:3); arrayfun(@(i) xSparseAll{c}(rezMax.dense.or(i), rezMax.dense.pos(i)), 1:3)];

    rezMax.sparse = struct();
    rezMax.sparse.bestBarStretch = stretchFactors(d);
    rezMax.sparse.bestLength = round(lenBarSparse * stretchFactors(d));
    [rezMax.sparse.maxcoef, rezMax.sparse.pos, rezMax.sparse.or] = get_best_parameters(zSparseAll{d}, 3, rezMax.sparse.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.sparse.maxcoefParts = [arrayfun(@(i) zDenseAll{d}(rezMax.sparse.or(i), rezMax.sparse.pos(i)), 1:3); arrayfun(@(i) zSparseAll{d}(rezMax.sparse.or(i), rezMax.sparse.pos(i)), 1:3)];
    rezMax.sparse.maxcoefPartsCC = [arrayfun(@(i) xDenseAll{d}(rezMax.sparse.or(i), rezMax.sparse.pos(i)), 1:3); arrayfun(@(i) xSparseAll{d}(rezMax.sparse.or(i), rezMax.sparse.pos(i)), 1:3)];

    if doExternal
      [~, e] = ismember(1, stretchFactors);
      rezMax.external = struct();
      rezMax.external.pos = externalAlignmentStruct.pos(idx);
      rezMax.external.or = externalAlignmentStruct.or(idx);
      rezMax.external.maxcoef = nansum([zDenseAll{e}(rezMax.external.or, rezMax.external.pos), zSparseAll{e}(rezMax.external.or, rezMax.external.pos)]);
      rezMax.external.maxcoefParts = [zDenseAll{e}(rezMax.external.or, rezMax.external.pos); zSparseAll{e}(rezMax.external.or, rezMax.external.pos)];
      rezMax.external.maxcoefPartsCC = [xDenseAll{e}(rezMax.external.or, rezMax.external.pos); xSparseAll{e}(rezMax.external.or, rezMax.external.pos)];
      rezMax.external.bestBarStretch = 1;
      rezMax.external.bestLength = lenBarSparse;
      fprintf("Barcode:\t%.0f\n\tS-Combined:\t%.2f\n\tS-Dense:\t%.2f\n\tS-Sparse:\t%.2f\n\tS-External:\t%.2f\n", idx, rezMax.dual.maxcoef(1), rezMax.dense.maxcoef(1), rezMax.sparse.maxcoef(1), rezMax.external.maxcoef)
    else
      fprintf("Barcode:\t%.0f\n\tS-Combined:\t%.2f\n\tS-Dense:\t%.2f\n\tS-Sparse:\t%.2f\n", idx, rezMax.dual.maxcoef(1), rezMax.dense.maxcoef(1), rezMax.sparse.maxcoef(1))
    end
    
    fprintf("\tPos. Offset:\t%.0f\n\tStr. Offset:\t%.0f\n", rezMax.dual.posOffset(1), rezMax.dual.stretchOffset(1))

    rezMaxM{idx} = rezMax;

  end
