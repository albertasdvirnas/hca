function rezMaxM = on_compare_dual_label(barcodeGen, theoryStruct, sets, externalAlignmentStruct)
  import CBT.Hca.Import.load_pval_struct;
  import CBT.Hca.UI.Helper.get_best_parameters;
  import DL.Hca.compute_zval
  
  if isfield(externalAlignmentStruct, 'molIds')
    doExternal = true;
  else
    doExternal = false;
  end

  %% Temp hardcoded settings
  digits(64) % Number of digits precision.

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
  psfPxRounded = round(sets.theory.psfSigmaWidth_nm / sets.theory.pixelWidth_nm, 4);

  try
    [ch1PsfInd, ch1Params] = load_pval_struct(sets.duallabel.pathZeroModelCh1);
    ch1ParamsInList = ismember(psfPxRounded, ch1PsfInd);
    ch1Params = ch1Params{ch1ParamsInList};
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", psfPxRounded, ch1ParamsPath));
  end

  try
    [ch2PsfInd, ch2Params] = load_pval_struct(sets.duallabel.pathZeroModelCh2);
    ch2ParamsInList = ismember(psfPxRounded, ch2PsfInd);
    ch2Params = ch2Params{ch2ParamsInList};
  catch
    error(compose("Zero-model params for psf: %.2f not found in file: %s", psfPxRounded, ch2ParamsPath));
  end

  %% Calculate scores
  numBarcodes = length(barcodeGen);
  rezMaxM = cell(1, numBarcodes);
  stretchFactors = sets.theory.stretchFactors;
  numStretchFactors = length(stretchFactors);
  numPxMatchDiff = sets.duallabel.numPxMatchDiff;
  stretchTolerance = sets.duallabel.stretchTolerance;
  positionTolerance = sets.duallabel.positionTolerance;
  barcodeGenCh1 = barcodeGen(1, :);
  barcodeGenCh2 = barcodeGen(2, :);

  for idx = 1:numBarcodes

    rezMax = struct();

    xCh1All = cell(1, numStretchFactors);
    xCh2All = cell(1, numStretchFactors);
    zCh1All = cell(1, numStretchFactors);
    zCh2All = cell(1, numStretchFactors);

    ch1Bar = barcodeGenCh1{idx}.rawBarcode;
    ch2Bar = barcodeGenCh2{idx}.rawBarcode;
    ch1Bit = barcodeGenCh1{idx}.rawBitmask;
    ch2Bit = barcodeGenCh2{idx}.rawBitmask;
    lenBarCh1 = length(ch1Bar);
    lenBarCh2 = length(ch2Bar);

    % run the loop for the stretch factors
    for j = 1:numStretchFactors
      
      % here interpolate both barcode and bitmask
      lenBarCh1Stretched = round(lenBarCh1 * stretchFactors(j));
      lenBarCh2Stretched = round(lenBarCh2 * stretchFactors(j));
      vCh1 = linspace(1, lenBarCh1, lenBarCh1Stretched);
      vCh2 = linspace(1, lenBarCh2, lenBarCh2Stretched);
      thisBarCh1 = interp1(ch1Bar, vCh1);
      thisBarCh2 = interp1(ch2Bar, vCh2);
      thisBitCh1 = ch1Bit(round(vCh1));
      thisBitCh2 = ch2Bit(round(vCh2));

      % Ch1 scores
      [zCh1, xCh1] = compute_zval(thisBarCh1, thisBitCh1, theorBar1, theorBitmask, sets.duallabel.barTypeCh1, ch1Params, sets.duallabel.ccThreshCh1, sets);

      % Ch2 scores
      [zCh2, xCh2] = compute_zval(thisBarCh2, thisBitCh2, theorBar2, theorBitmask, sets.duallabel.barTypeCh2, ch2Params, sets.duallabel.ccThreshCh2, sets);

      xCh1All{j} = xCh1;
      xCh2All{j} = xCh2;
      zCh1All{j} = zCh1;
      zCh2All{j} = zCh2;

    end

    % Optimise combined score
    zDualAll = cell(1, numStretchFactors);
    posOffset = cell(1, numStretchFactors);
    stretchOffset = cell(1, numStretchFactors);

    for j = 1:numStretchFactors
      stretchIds = max(1, j - stretchTolerance):min(numStretchFactors, j + stretchTolerance);

      thisCh2Z = zCh2All{j};
      thisDualZ = zeros([size(zCh2All{j}), length(stretchIds), 2 * positionTolerance + 1]);

      for k = 1:length(stretchIds)
        for m = 1:2 * positionTolerance + 1
          shift = -(m - positionTolerance - 1);
          thisDualZ(:, :, k, m) = nansum(cat(3, thisCh2Z, circshift(zCh1All{stretchIds(k)}, [0 shift])), 3);
        end
      end

      tmpScoreMat1 = nanmax(thisDualZ, [], 4);
      tmpScoreMat2 = nanmax(thisDualZ, [], 3);
      [~, tmpStretch] = nanmax(tmpScoreMat1, [], 3);
      [zDualAll{j}, tmpPos] = nanmax(tmpScoreMat2, [], 4);
      stretchOffset{j} = reshape(stretchIds(tmpStretch) - j, size(zCh2All{j}));
      posOffset{j} = reshape(tmpPos - positionTolerance - 1, size(zCh2All{j}));
      
      isUnchanged = zDualAll{j} == zCh1All{j} | zDualAll{j} == zCh2All{j};
      stretchOffset{j}(isUnchanged) = 0;
      posOffset{j}(isUnchanged) = 0;
    end

    % Get best stretchF for individual barcodes
    tmpLongestBar = round(max(lenBarCh1, lenBarCh2) * stretchFactors(end));
    [~, b] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zDualAll));
    [~, c] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zCh1All));
    [~, d] = nanmax(cellfun(@(x) get_best_parameters(x, 1, tmpLongestBar, 1, numPxMatchDiff, theorUnkCount), zCh2All));

    rezMax.dual = struct();
    rezMax.dual.bestBarStretch = stretchFactors(b);
    rezMax.dual.bestLength = round(max(lenBarCh1, lenBarCh2) * stretchFactors(b));
    [rezMax.dual.maxcoef, rezMax.dual.pos, rezMax.dual.or] = get_best_parameters(zDualAll{b}, 3, rezMax.dual.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.dual.posOffset = arrayfun(@(i) posOffset{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3);
    rezMax.dual.stretchOffset = arrayfun(@(i) stretchOffset{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3);
    optPos = max(1, rezMax.dual.pos + rezMax.dual.posOffset);
    optStr = b + rezMax.dual.stretchOffset;
    rezMax.dual.optStr = stretchFactors(optStr);
    rezMax.dual.maxcoefParts = [arrayfun(@(i) zCh1All{optStr(i)}(rezMax.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) zCh2All{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3)];
    rezMax.dual.maxcoefPartsCC = [arrayfun(@(i) xCh1All{optStr(i)}(rezMax.dual.or(i), optPos(i)), 1:3); arrayfun(@(i) xCh2All{b}(rezMax.dual.or(i), rezMax.dual.pos(i)), 1:3)];

    rezMax.ch1 = struct();
    rezMax.ch1.bestBarStretch = stretchFactors(c);
    rezMax.ch1.bestLength = round(lenBarCh1 * stretchFactors(c));
    [rezMax.ch1.maxcoef, rezMax.ch1.pos, rezMax.ch1.or] = get_best_parameters(zCh1All{c}, 3, rezMax.ch1.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.ch1.maxcoefParts = [arrayfun(@(i) zCh1All{c}(rezMax.ch1.or(i), rezMax.ch1.pos(i)), 1:3); arrayfun(@(i) zCh2All{c}(rezMax.ch1.or(i), rezMax.ch1.pos(i)), 1:3)];
    rezMax.ch1.maxcoefPartsCC = [arrayfun(@(i) xCh1All{c}(rezMax.ch1.or(i), rezMax.ch1.pos(i)), 1:3); arrayfun(@(i) xCh2All{c}(rezMax.ch1.or(i), rezMax.ch1.pos(i)), 1:3)];

    rezMax.ch2 = struct();
    rezMax.ch2.bestBarStretch = stretchFactors(d);
    rezMax.ch2.bestLength = round(lenBarCh2 * stretchFactors(d));
    [rezMax.ch2.maxcoef, rezMax.ch2.pos, rezMax.ch2.or] = get_best_parameters(zCh2All{d}, 3, rezMax.ch2.bestLength, 1, numPxMatchDiff, theorUnkCount);
    rezMax.ch2.maxcoefParts = [arrayfun(@(i) zCh1All{d}(rezMax.ch2.or(i), rezMax.ch2.pos(i)), 1:3); arrayfun(@(i) zCh2All{d}(rezMax.ch2.or(i), rezMax.ch2.pos(i)), 1:3)];
    rezMax.ch2.maxcoefPartsCC = [arrayfun(@(i) xCh1All{d}(rezMax.ch2.or(i), rezMax.ch2.pos(i)), 1:3); arrayfun(@(i) xCh2All{d}(rezMax.ch2.or(i), rezMax.ch2.pos(i)), 1:3)];

    if doExternal
      [~, e] = ismember(1, stretchFactors);
      rezMax.external = struct();
      rezMax.external.pos = externalAlignmentStruct.pos(idx);
      rezMax.external.or = externalAlignmentStruct.or(idx);
      rezMax.external.maxcoef = nansum([zCh1All{e}(rezMax.external.or, rezMax.external.pos), zCh2All{e}(rezMax.external.or, rezMax.external.pos)]);
      rezMax.external.maxcoefParts = [zCh1All{e}(rezMax.external.or, rezMax.external.pos); zCh2All{e}(rezMax.external.or, rezMax.external.pos)];
      rezMax.external.maxcoefPartsCC = [xCh1All{e}(rezMax.external.or, rezMax.external.pos); xCh2All{e}(rezMax.external.or, rezMax.external.pos)];
      rezMax.external.bestBarStretch = 1;
      rezMax.external.bestLength = lenBarCh2;
      fprintf("Barcode:\t%.0f\n\tS-Combined:\t%.2f\n\tS-Channel1:\t%.2f\n\tS-Channel2:\t%.2f\n\tS-External:\t%.2f\n", idx, rezMax.dual.maxcoef(1), rezMax.ch1.maxcoef(1), rezMax.ch2.maxcoef(1), rezMax.external.maxcoef)
    else
      fprintf("Barcode:\t%.0f\n\tS-Combined:\t%.2f\n\tS-Channel1:\t%.2f\n\tS-Channel2:\t%.2f\n", idx, rezMax.dual.maxcoef(1), rezMax.ch1.maxcoef(1), rezMax.ch2.maxcoef(1))
    end
    
    fprintf("\tPos. Offset:\t%.0f\n\tStr. Offset:\t%.0f\n", rezMax.dual.posOffset(1), rezMax.dual.stretchOffset(1))

    rezMaxM{idx} = rezMax;

  end
