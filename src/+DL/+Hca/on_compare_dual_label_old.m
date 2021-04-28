function rezMaxM = on_compare_dual_label(...
    barcodeGen, ...
    theoryStruct, ...
    numPixelsAroundBestTheoryMask, ...
    sets, ...
    doBitmask, ...
    externalAlignmentStruct)

  % Change to balance precision/speed
  digits(128) % Number of digits precision.
  xThreshDense = 0; % Minimum value of CC to calculate score.
  xThreshSparse = 0; % Minimum value of CC to calculate score.
  zThresh = -inf;
  %   stretchTol = 5;
  %   posTol = 10;

  % load theory barcode txt file.
  formatSpec = '%f';
  fileID = fopen(theoryStruct.filename, 'r');
  theorBar1 = transpose(fscanf(fileID, formatSpec));
  fclose(fileID);
  fileID = fopen(theoryStruct.filename2, 'r');
  theorBar2 = transpose(fscanf(fileID, formatSpec));
  fclose(fileID);

  try
    fileID = fopen(strrep(theoryStruct.filename, 'barcode', 'bitmask'), 'r');
    theorBit = transpose(fscanf(fileID, formatSpec));
    fclose(fileID);
  catch
    theorBit = ones(size(theorBar1));
  end

  % Load zero-model params
  theoryNumPixels = 3088286401 / sets.pvalue.pixelWidth_nm * sets.pvalue.nmbp;
  numStretchFactors = length(sets.theory.stretchFactors);
  import CBT.Hca.Import.load_pval_struct;
  cbParamsPath = fullfile(sets.duallabel.paramFolder, [sets.duallabel.denseBarcodeType "_zero_model_params.txt"]);
  dotsParamsPath = fullfile(sets.duallabel.paramFolder, compose("dots_%s_zero_model_params.txt", sets.theory.sparsePattern));
  [cbPsfInd, cbParams] = load_pval_struct(cbParamsPath);
  [dotsPsfInd, dotParams] = load_pval_struct(dotsParamsPath);
  psfPxRounded = round(sets.pvalue.psfSigmaWidth_nm / sets.pvalue.pixelWidth_nm, 4);

  try
    cbParamLambda = cbParams{ismember(psfPxRounded, cbPsfInd)}(1);
    cbParamExa = cbParams{ismember(psfPxRounded, cbPsfInd)}(2);
    cbParamNu = cbParams{ismember(psfPxRounded, cbPsfInd)}(3);
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", ...
      psfPxRounded, cbParamsPath));
  end

  try
    dotParamLambda = dotParams{ismember(psfPxRounded, dotsPsfInd)}(1);
    dotParamExa = dotParams{ismember(psfPxRounded, dotsPsfInd)}(2);
    dotParamMu = dotParams{ismember(psfPxRounded, dotsPsfInd)}(3);
    dotParamXi = dotParams{ismember(psfPxRounded, dotsPsfInd)}(4);
    dotParamSigma = dotParams{ismember(psfPxRounded, dotsPsfInd)}(5);
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", ...
      psfPxRounded, dotsParamsPath));
  end

  rezMaxM = cell(1, size(barcodeGen, 2));
  barcodeGenDense = barcodeGen(1, :);
  barcodeGenSparse = barcodeGen(2, :);

  % for all the barcodes run
  for idx = 1:size(barcodeGen, 2)
    tic

    % rezMaz stores the results for one barcode
    %     xDenseAll = cell(1, numStretchFactors);
    %     xSparseAll = cell(1, numStretchFactors);
    %     zDenseAll = cell(1, numStretchFactors);
    %     zSparseAll = cell(1, numStretchFactors);
    rezMaxDual = cell(1, numStretchFactors);
    rezMaxDense = cell(1, numStretchFactors);
    rezMaxSparse = cell(1, numStretchFactors);
    rezMaxExternal = cell(1, numStretchFactors);

    % length of this barcode
    lenBarTested = length(barcodeGenDense{idx}.rawBarcode);

    bitmask = barcodeGenDense{idx}.rawBitmask;

    % run the loop for the stretch factors
    for j = 1:numStretchFactors
      % here interpolate both barcode and bitmask
      lenBarStretched = round(lenBarTested * sets.theory.stretchFactors(j));
      v = linspace(1, lenBarTested, lenBarStretched);
      barDense = interp1(barcodeGenDense{idx}.rawBarcode, v);
      barSparse = interp1(barcodeGenSparse{idx}.rawBarcode, v);
      thisBitmask = bitmask(round(v));

      import DL.Hca.beta_ev_cdf

      if nargin > 4 && doBitmask
        nuLenBar = sum(thisBitmask);
        import SignalRegistration.masked_pcc_corr
        xDense = masked_pcc_corr(barDense, theorBar1, thisBitmask, true(size(theorBar1)));
      else
        xDense = MASS_PCC(...
          theorBar1, ...
          barDense, ...
          2^(4 + nextpow2(lenBarTested)));
        nuLenBar = lenBarStretched;
      end

      zDense = zThresh * ones(size(xDense));
      lambdaEffDense = cbParamLambda * cbParamExa * numStretchFactors * (theoryNumPixels - lenBarStretched);
      doZScore = xDense > xThreshDense;
      zDense(doZScore) = -sqrt(2) * erfcinv(2 * beta_ev_cdf(...
        xDense(doZScore), ...
        max(4, cbParamNu * nuLenBar), ...
        1, ...
        lambdaEffDense));
      doIncreasedPrecision = isinf(zDense) & zDense > 0;
      zDense(doIncreasedPrecision) = -sqrt(2) * erfcinv(2 * beta_ev_cdf(...
        xDense(doIncreasedPrecision), ...
        max(4, cbParamNu * nuLenBar), ...
        1, ...
        lambdaEffDense, ...
        true));

      import DL.Hca.trunc_normal_ev_cdf
      xSparse = MASS_DOT_CC(...
        theorBar2, ...
        barSparse, ...
        2^(4 + nextpow2(lenBarTested)));
      zSparse = zThresh * ones(size(xSparse));
      lambdaEffSparse = dotParamLambda * dotParamExa * numStretchFactors * (theoryNumPixels - lenBarStretched);
      doZScore = xSparse > xThreshSparse;
      zSparse(doZScore) = -sqrt(2) * erfcinv(2 * trunc_normal_ev_cdf(...
        xSparse(doZScore), ...
        -sqrt(2) * erfcinv(2 * linspace(0.001, 0.999, 31), ...
        dotParamMu, ...
        1 / sqrt(dotParamXi * lenBarStretched)), ...
        1 / sqrt(dotParamSigma * lenBarStretched), ...
        lambdaEffSparse, ...
        0, 1));
      doIncreasedPrecision = isinf(zSparse) & zSparse > 0;
      zSparse(doIncreasedPrecision) = -sqrt(2) * erfcinv(2 * trunc_normal_ev_cdf(...
        xSparse(doIncreasedPrecision), ...
        -sqrt(2) * erfcinv(2 * linspace(0.001, 0.999, 51), ...
        dotParamMu, ...
        1 / sqrt(dotParamXi * lenBarStretched)), ...
        1 / sqrt(dotParamSigma * lenBarStretched), ...
        lambdaEffSparse, ...
        0, 1, true));

      %       xDenseAll{j} = xDense;
      %       xSparseAll{j} = xSparse;
      %       zDenseAll{j} = zDense;
      %       zSparseAll{j} = zSparse;

      import CBT.Hca.UI.Helper.get_best_parameters;
      [rezMaxDual{j}.maxcoef, rezMaxDual{j}.pos, rezMaxDual{j}.or] = get_best_parameters(...
        (zDense + zSparse) / sqrt(2), 3, lenBarStretched, 1, numPixelsAroundBestTheoryMask, theorBit);
      [rezMaxDense{j}.maxcoef, rezMaxDense{j}.pos, rezMaxDense{j}.or] = get_best_parameters(...
        zDense, 3, lenBarStretched, 1, numPixelsAroundBestTheoryMask, theorBit);
      [rezMaxSparse{j}.maxcoef, rezMaxSparse{j}.pos, rezMaxSparse{j}.or] = get_best_parameters(...
        zSparse, 3, lenBarStretched, 1, numPixelsAroundBestTheoryMask, theorBit);

      rezMaxDual{j}.maxcoefParts = [arrayfun(@(i) zDense(rezMaxDual{j}.or(i), rezMaxDual{j}.pos(i)), 1:3); arrayfun(@(i) zSparse(rezMaxDual{j}.or(i), rezMaxDual{j}.pos(i)), 1:3)];
      rezMaxDual{j}.maxcoefPartsCC = [arrayfun(@(i) xDense(rezMaxDual{j}.or(i), rezMaxDual{j}.pos(i)), 1:3); arrayfun(@(i) xSparse(rezMaxDual{j}.or(i), rezMaxDual{j}.pos(i)), 1:3)];
      rezMaxDense{j}.maxcoefParts = [arrayfun(@(i) zDense(rezMaxDense{j}.or(i), rezMaxDense{j}.pos(i)), 1:3); arrayfun(@(i) zSparse(rezMaxDense{j}.or(i), rezMaxDense{j}.pos(i)), 1:3)];
      rezMaxDense{j}.maxcoefPartsCC = [arrayfun(@(i) xDense(rezMaxDense{j}.or(i), rezMaxDense{j}.pos(i)), 1:3); arrayfun(@(i) xSparse(rezMaxDense{j}.or(i), rezMaxDense{j}.pos(i)), 1:3)];
      rezMaxSparse{j}.maxcoefParts = [arrayfun(@(i) zDense(rezMaxSparse{j}.or(i), rezMaxSparse{j}.pos(i)), 1:3); arrayfun(@(i) zSparse(rezMaxSparse{j}.or(i), rezMaxSparse{j}.pos(i)), 1:3)];
      rezMaxSparse{j}.maxcoefPartsCC = [arrayfun(@(i) xDense(rezMaxSparse{j}.or(i), rezMaxSparse{j}.pos(i)), 1:3); arrayfun(@(i) xSparse(rezMaxSparse{j}.or(i), rezMaxSparse{j}.pos(i)), 1:3)];

      if nargin > 5
        rezMaxExternal{j}.pos = externalAlignmentStruct.pos(idx);
        rezMaxExternal{j}.or = externalAlignmentStruct.or(idx);
        rezMaxExternal{j}.maxcoef = (zDense(rezMaxExternal{j}.or, rezMaxExternal{j}.pos) + zSparse(rezMaxExternal{j}.or, rezMaxExternal{j}.pos)) / sqrt(2);
        rezMaxExternal{j}.maxcoefParts = [zDense(rezMaxExternal{j}.or, rezMaxExternal{j}.pos); zSparse(rezMaxExternal{j}.or, rezMaxExternal{j}.pos)];
        rezMaxExternal{j}.maxcoefPartsCC = [xDense(rezMaxExternal{j}.or, rezMaxExternal{j}.pos); xSparse(rezMaxExternal{j}.or, rezMaxExternal{j}.pos)];
      else
        rezMaxExternal{j}.maxcoef = nan;
      end

    end

    % Optimise combined score
    %     zDual = cell(1, numStretchFactors);
    %     posOffset = cell(1, numStretchFactors);
    %     strOffset = cell(1, numStretchFactors);
    %
    %     for j = 1:numStretchFactors
    %       tmpScoreMat2 = zThresh*ones([size(zSparseAll{j}) 2 * stretchTol + 1]);
    %       tmpPosOffset = (posTol + 1)*ones(size(tmpScoreMat2));
    %
    %       for n = 1:2 * stretchTol + 1
    %         stretchId = j + n - stretchTol - 1;
    %
    %         if stretchId < 1 || stretchId > numStretchFactors; continue; end
    %
    %         tmpScoreMat = zThresh*ones([size(zSparseAll{j}) 2 * posTol + 1]);
    %         for m = 1:2 * posTol + 1
    %           shift = -(m - posTol - 1);
    %           tmpScoreMat(:, :, m) = (zSparseAll{j} + circshift(zDenseAll{stretchId}, [0 shift])) / sqrt(2);
    %         end
    %
    %         [tmpScoreMat2(:, :, n), tmpPosOffset(:, :, n)] = nanmax(tmpScoreMat, [], 3);
    %       end
    %
    %       [zDual{j}, strOffset{j}] = nanmax(tmpScoreMat2, [], 3);
    %       [tmprow, tmpcol] = ind2sub(size(zSparseAll{j}), 1:numel(zSparseAll{j}));
    %       posOffset{j} = reshape(tmpPosOffset(sub2ind(size(tmpPosOffset), tmprow, tmpcol, strOffset{j}(:)')), size(zSparseAll{j}));
    %     end

    % Get best stretchF for individual barcodes
    %     [~, b] = nanmax(cellfun(@(x) nanmax(x, [], 'all'), zDual));
    [~, b] = nanmax(cellfun(@(x) x.maxcoef(1), rezMaxDual));
    [~, c] = nanmax(cellfun(@(x) x.maxcoef(1), rezMaxDense));
    [~, d] = nanmax(cellfun(@(x) x.maxcoef(1), rezMaxSparse));
    [~, e] = nanmax(cellfun(@(x) x.maxcoef(1), rezMaxExternal));

    %     rezMaxM{idx}.dual = struct();
    %     rezMaxM{idx}.dual.bestBarStretch = sets.theory.stretchFactors(b);
    %     rezMaxM{idx}.dual.bestLength = round(lenBarTested * sets.theory.stretchFactors(b));
    %     import CBT.Hca.UI.Helper.get_best_parameters;
    %     [rezMaxM{idx}.dual.maxcoef, rezMaxM{idx}.dual.pos, rezMaxM{idx}.dual.or] = get_best_parameters(...
    %       zDual{b}, 3, rezMaxM{idx}.dual.bestLength, 1, numPixelsAroundBestTheoryMask, theorBit);
    %     optPosOffset = arrayfun(@(i) posOffset{b}(rezMaxM{idx}.dual.or(i), rezMaxM{idx}.dual.pos(i)), 1:3);
    %     optStrOffset = arrayfun(@(i) strOffset{b}(rezMaxM{idx}.dual.or(i), rezMaxM{idx}.dual.pos(i)), 1:3);
    %     optPos = max(1, rezMaxM{idx}.dual.pos + optPosOffset - posTol - 1);
    %     optStr = max(1, b + optStrOffset - stretchTol - 1);
    %     rezMaxM{idx}.dual.maxcoefParts = [arrayfun(@(i) zDenseAll{optStr(i)}(rezMaxM{idx}.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) zSparseAll{b}(rezMaxM{idx}.dual.or(i), rezMaxM{idx}.dual.pos(i)), 1:3)];
    %     rezMaxM{idx}.dual.maxcoefPartsCC = [arrayfun(@(i) xDenseAll{optStr(i)}(rezMaxM{idx}.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) xSparseAll{b}(rezMaxM{idx}.dual.or(i), rezMaxM{idx}.dual.pos(i)), 1:3)];
    rezMaxM{idx}.dual = rezMaxDual{b};
    rezMaxM{idx}.dual.bestBarStretch = sets.theory.stretchFactors(b);
    rezMaxM{idx}.dual.bestLength = round(lenBarTested * sets.theory.stretchFactors(b));
    rezMaxM{idx}.dense = rezMaxDense{c};
    rezMaxM{idx}.dense.bestBarStretch = sets.theory.stretchFactors(c);
    rezMaxM{idx}.dense.bestLength = round(lenBarTested * sets.theory.stretchFactors(c));
    rezMaxM{idx}.sparse = rezMaxSparse{d};
    rezMaxM{idx}.sparse.bestBarStretch = sets.theory.stretchFactors(d);
    rezMaxM{idx}.sparse.bestLength = round(lenBarTested * sets.theory.stretchFactors(d));
    rezMaxM{idx}.external = rezMaxExternal{e};
    rezMaxM{idx}.external.bestBarStretch = sets.theory.stretchFactors(e);
    rezMaxM{idx}.external.bestLength = round(lenBarTested * sets.theory.stretchFactors(e));

    if nargin > 5
      disp(num2str([idx rezMaxM{idx}.dual.maxcoef(1) rezMaxM{idx}.dense.maxcoef(1) rezMaxM{idx}.sparse.maxcoef(1) rezMaxM{idx}.external.maxcoef rezMaxM{idx}.external.maxcoefParts']))
    else
      disp(num2str([idx rezMaxM{idx}.dual.maxcoef(1) rezMaxM{idx}.dense.maxcoef(1) rezMaxM{idx}.sparse.maxcoef(1)]))
    end

    toc
  end
