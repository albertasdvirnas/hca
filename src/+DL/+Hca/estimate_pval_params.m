function params = estimate_pval_params(sets)
%% estimate_pval_params
import DL.Hca.beta_ev_fit
import DL.Hca.trunc_normal_ev_fit

% Load settings
lenLong = sets.pvalue.lenLong;
lenShort = sets.pvalue.lenShort;
psfPx = sets.pvalue.psfSigmaWidth_nm/sets.pvalue.pixelWidth_nm;
stretchMax = 0; % No reason to collect scores for any other value.
stretchStep = 0.01; % Temp?
lenLong_name = strcat(num2str(round(lenLong/10^floor(log10(lenLong))), 1), ...
  "e", num2str(floor(log10(lenLong))));
lenShort_name = num2str(lenShort);
psfPx_name = num2str(psfPx, 3);
stretchMax_name = num2str(stretchMax);

% Expected file names and paths
longVarFileName = compose("%s_scores_long_%s_short_%s_psf_%s_stretch_%s.txt", ...
  sets.pvalue.barcodeType, "var", lenShort_name, psfPx_name, stretchMax_name);
longVarFilePath = fullfile(sets.pvalue.scoreFolder, longVarFileName);
shortVarFileName = compose("%s_scores_long_%s_short_%s_psf_%s_stretch_%s.txt", ...
  sets.pvalue.barcodeType, lenLong_name, "var", psfPx_name, stretchMax_name);
shortVarFilePath = fullfile(sets.pvalue.scoreFolder, shortVarFileName);
stretchVarFileName = compose("%s_scores_long_%s_short_%s_psf_%s_stretch_%s.txt", ...
  sets.pvalue.barcodeType, lenLong_name, lenShort_name, psfPx_name, "var");
stretchVarFilePath = fullfile(sets.pvalue.scoreFolder, stretchVarFileName);
import CBT.Hca.Import.load_pval_struct;
[longVar, longVarData] = load_pval_struct(longVarFilePath);
if isempty(longVar)
  error("File containing scores with variable long barcode length not found.")
end
[shortVar, shortVarData] = load_pval_struct(shortVarFilePath);
if isempty(shortVar)
  error("File containing scores with variable short barcode length not found.")
end
[stretchVar, stretchVarData] = load_pval_struct(stretchVarFilePath);
if isempty(stretchVar)
  error("File containing scores with variable stretch factor not found.")
end
if not(isfolder(sets.pvalue.paramFolder))
  mkdir(sets.pvalue.paramFolder);
end
switch sets.pvalue.barcodeType
  case 'cb'
    paramFileName = "cb_zero_model_params.txt";
  case 'dots'
    paramFileName = compose("dots_%s_zero_model_params.txt", sets.pvalue.pattern);
end
paramFilePath = fullfile(sets.pvalue.paramFolder, paramFileName);
import CBT.Hca.Import.load_pval_struct;
[vals, data] = load_pval_struct(paramFilePath);
if ismember(sets.pvalue.psfSigmaWidth_nm, vals)
  error(compose("Zero-model parameters for this PSF are already in file: %s", paramFilePath))
end

%% Estimate Lambda scaling by long barcode length.
% Start by finding the mean of the params that are expected to be constant,
% then use them to find the linear fit of Lambda.
disp("Estimating Lambda scaling by long barcode length")
longVarLambda = nan(length(longVar), 1);
switch sets.pvalue.barcodeType
  case 'cb'
    longVarNu = nan(length(longVar), 1);
    for i=1:length(longVar)
      longVarNu(i) = beta_ev_fit( ...
        longVarData{i}, ...
        [4 1 1], ...
        [inf 1 2*(longVar(i)-lenShort)], ...
        [4 1 1], ...
        [false true false]);
    end
    for i=1:length(longVar)
      [~, ~, longVarLambda(i)] = beta_ev_fit( ...
        longVarData{i}, ...
        [4 1 1], ...
        [inf 1 2*(longVar(i)-lenShort)], ...
        [nanmean(longVarNu) 1 1], ...
        [true true false]);
    end
  case 'dots'
    % TODO
end
fitLambda = fit( ...
  longVar(:), ...
  longVarLambda, ...
  'poly1', ...
  'Lower', [0 0], ...
  'Upper', [inf 0]);
params(1) = fitLambda.p1;

%% Estimate Lambda scaling by stretch factor
% Uses some of the results from the previous Lambda fit.
disp("Estimating Lambda scaling by stretch factor")
stretchVarExa = nan(length(stretchVar), 1);
switch sets.pvalue.barcodeType
  case 'cb'
    for i=1:length(stretchVar)
      [~, ~, stretchVarExa(i)] = beta_ev_fit( ...
        stretchVarData{i}, ...
        [4 1 1], ...
        [inf 1 inf], ...
        [nanmean(longVarNu) 1 1], ...
        [true true false]);
    end
  case 'dots'
    % TODO
end
fitExa = fit( ...
  2*stretchVar(:)/stretchStep+1, ...
  stretchVarExa/(params(1)*(lenLong-lenShort)), ...
  'poly1', ...
  'Lower', [0 0], ...
  'Upper', [inf 0]);
params(2) = fitExa.p1;

%% Estimate distribution specific params scaling by short barcode length.
% Uses some of the results from the first Lambda fit.
disp("Estimating distribution specific params scaling by short barcode length")
switch sets.pvalue.barcodeType
  case 'cb'
    shortVarNu = nan(length(shortVar), 1);
    for i=1:length(shortVar)
      shortVarNu(i) = beta_ev_fit( ...
        shortVarData{i}, ...
        [4 1 1], ...
        [inf 1 2*(lenLong-shortVar(i))], ...
        [nanmean(longVarNu)/lenShort*shortVar(i) 1 params(1)*(lenLong-shortVar(i))], ...
        [false true true]);
    end
    fitNu = fit( ...
      shortVar(:), ...
      shortVarNu, ...
      'poly1', ...
      'Lower', [0 0], ...
      'Upper', [inf 0]);
    params(3) = fitNu.p1;
  case 'dots'
    % TODO
end

%% Write output to file
import CBT.Hca.Export.export_pval_struct;
vals(end+1) = sets.pvalue.psfSigmaWidth_nm;
data{end+1} = params;
export_pval_struct(paramFilePath, vals, data);
fprintf("Finished estimating params, wrote output to: %s\n", paramFilePath)