function [data, means] = run_random_score_gen(sets)
%% run_random_score_gen

rng('shuffle')

% Load settings
bppx = sets.pvalue.pixelWidth_nm/sets.pvalue.nmbp;
lenLong_all = sets.pvalue.lenLong;
lenShort_all = sets.pvalue.lenShort;
psfPx_all = sets.pvalue.psfSigmaWidth_nm/sets.pvalue.pixelWidth_nm;
stretchMax_all = sets.pvalue.stretchMax;
lenLong_name = strcat(num2str(round(lenLong_all/10^floor(log10(lenLong_all))), 1), ...
  "e", num2str(floor(log10(lenLong_all))));
lenShort_name = num2str(lenShort_all);
psfPx_name = num2str(psfPx_all, 3);
stretchMax_name = num2str(stretchMax_all);
switch sets.pvalue.variableType
  case 'lenLong'
    lenLong_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    lenLong_name = "var";
    vals = lenLong_all;
  case 'lenShort'
    lenShort_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    lenShort_name = "var";
    vals = lenShort_all;
  case 'psf'
    psfPx_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    psfPx_name = "var";
    vals = psfPx_all;
  case 'stretchMax'
    stretchMax_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    stretchMax_name = "var";
    vals = stretchMax_all;
  otherwise
    error('Bad variable type selected for random score generation')
end

% Input/Output
if not(isfolder(sets.pvalue.scoreFolder))
  mkdir(sets.pvalue.scoreFolder);
end
if not(isfolder(sets.pvalue.fastaFolder))
  mkdir(sets.pvalue.fastaFolder);
end
if strcmp(sets.pvalue.barcodeType, 'dots')
  typeName = sets.pvalue.pattern;
else
  typeName = sets.pvalue.barcodeType;
end
scoreFileName = compose("%s_scores_long_%s_short_%s_psf_%s_stretch_%s.txt", ...
  typeName, lenLong_name, lenShort_name, psfPx_name, stretchMax_name);
meansFileName = compose("%s_means_long_%s_short_%s_psf_%s_stretch_%s.txt", ...
  typeName, lenLong_name, lenShort_name, psfPx_name, stretchMax_name);
scoreFilePath = fullfile(sets.pvalue.scoreFolder, scoreFileName);
meansFilePath = fullfile(sets.pvalue.scoreFolder, meansFileName);
import CBT.Hca.Import.load_pval_struct;
[~, data] = load_pval_struct(scoreFilePath);
[~, means] = load_pval_struct(meansFilePath);
alreadyInData = length(data);
if strcmp(sets.pvalue.barcodeType, 'cb')
    switch sets.pvalue.cbRandomizationMethod, ...
      case {'phaseRandomization', 'hybrid'}
        [meanfftFilename, dirpath] = uigetfile('*.mat', 'Select file with mean fft magnitudes');
        if isequal(dirpath, 0)
            return;
        end
        meanfftFile = load(fullfile(dirpath, meanfftFilename));
        meanfft = meanfftFile.meanFFT;
      case 'sample'
        fid = fopen(sets.theories);
        thrNames = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        for i=1:length(thrNames{1})
          fid = fopen(thrNames{1}{i});
          tmpStream = textscan(fid,'%f','delimiter','\t');
          theories{i} = tmpStream{1}';
          fclose(fid);
        end
        theoryNumPixels = round(3088286401/sets.pvalue.pixelWidth_nm*sets.pvalue.nmbp);
    end
end

%%
import CBT.Hca.Core.Pvalue.convolve_bar;

% For each psf
for psfPx=psfPx_all
  psfBp = round(psfPx*bppx);
  % For each long length
  for lenLong=lenLong_all
    if strcmp(sets.pvalue.variableType, 'lenLong') && alreadyInData > 0
        alreadyInData = alreadyInData - 1;
        continue
    end
    % Generate long barcode
    lenLongBp = round(lenLong*bppx);
    switch sets.pvalue.barcodeType
      case 'cb'
        switch sets.pvalue.cbRandomizationMethod
          case 'randn'
            randLong = normrnd(0, 1, 1, lenLong);
            randLong = convolve_bar(randLong, psfPx, lenLong);
          case {'phaseRandomization', 'hybrid'}
            import CBT.RandBarcodeGen.PhaseRandomization.generate_rand_barcodes_from_fft_zero_model
            randLong = generate_rand_barcodes_from_fft_zero_model(meanfft, ...
              ceil(lenLong/length(meanfft)), length(meanfft));
            randLong = [randLong{:}];
            t = randi(length(randLong)-lenLong+1);
            randLong = randLong(t:t+lenLong-1);
            if strcmp(sets.pvalue.cbRandomizationMethod, 'hybrid')
              randLongExtra = normrnd(0, 1, 1, lenLong);
              randLong = 0.7*zscore(randLong) +  0.3*zscore(convolve_bar(randLongExtra, psfPx, lenLong));
            end
          case 'sample'
            flip = logical(randi(2, 1, length(theories))-1);
            theories(flip) = cellfun(@fliplr, theories(flip), 'un', 0);
            theories = theories(randperm(length(theories)));
            t = randi(theoryNumPixels-lenLong+1);
            cl = cumsum(cellfun(@length, theories));
            t1 = find(cl >= t, 1);
            t2 = find(cl >= t+lenLong-1, 1);
            randLong = [theories{t1:t2}];
            t = randi(length(randLong)-lenLong+1);
            randLong = randLong(t:t+lenLong-1);
            unusedTheories = [theories{1:t1-1} theories{t2+1:end}];
        end
      case {'dots', 'ecodam'}
        switch sets.pvalue.barcodeType
          case 'dots'
            pattern = sets.pvalue.pattern;
          case 'ecodam'
            pattern = 'GATC';
        end
        randFastaName = compose("random_fasta_length_%.0f.fasta", lenLongBp);
        randFastaPath = fullfile(sets.pvalue.fastaFolder, randFastaName);
        if not(isfile(randFastaPath))
          import DL.Hca.generate_random_fasta
          generate_random_fasta(randFastaPath, lenLongBp)
        end
        import DL.Hca.compute_theory_dotbar
        randLong = compute_theory_dotbar( ...
          randFastaPath, ...
          pattern, ...
          bppx, ...
          psfBp, ...
          2^15, 2^13);
        randLong = randLong(:)';
    end
    % For each stretch-max
    for stretchMax=stretchMax_all
      % For each short length
      for lenShort=lenShort_all
        % Skip data computed previously
        if alreadyInData > 0
          alreadyInData = alreadyInData - 1;
          continue
        end
        fprintf("=====================================================\n");
        fprintf("Generating random scores for params:\n")
        fprintf("Long length:\t%.0f\nShort Length:\t%.0f\nPsf width:\t\t%.2f\nStretch max:\t%.2f\n", lenLong, lenShort, psfPx, stretchMax);
        bestScoreVec = nan(1, sets.pvalue.numRnd);
        scoreMeanVec = nan(1, sets.pvalue.numRnd);
        f = waitbar( ...
          0, ...
          "Comparing barcode pair 1, estimated time remaining: inf minutes.", ...
            'Position', [0, 0, 200, 50]);
        tic
        % For each random sample
        for i=1:sets.pvalue.numRnd
          % Generate short barcode
          lenShortBp = round(lenShort*bppx);
          switch sets.pvalue.barcodeType
            case 'cb'
              switch sets.pvalue.cbRandomizationMethod
                case 'randn'
                  randShort = normrnd(0, 1, 1, lenShort);
                  randShort = convolve_bar(randShort, psfPx, lenShort);
                case {'phaseRandomization', 'hybrid'}
                  import CBT.RandBarcodeGen.PhaseRandomization.generate_rand_barcodes_from_fft_zero_model
                  randShort = generate_rand_barcodes_from_fft_zero_model(meanfft, ...
                    ceil(lenShort/length(meanfft)), length(meanfft));
                  randShort = [randShort{:}];
                  t = randi(length(randShort)-lenShort+1);
                  randShort = randShort(t:t+lenShort-1);
                  if strcmp(sets.pvalue.cbRandomizationMethod, 'hybrid')
                    randShortExtra = normrnd(0, 1, 1, lenShort);
                    randShort = 0.7*zscore(randShort) + 0.3*zscore(convolve_bar(randShortExtra, psfPx, lenShort));
                  end
                case 'sample'
                  t = randi(length(unusedTheories)-lenShort+1);
                  randShort = unusedTheories(t:t+lenShort-1);
              end
            case {'dots', 'ecodam'}
              switch sets.pvalue.barcodeType
                case 'dots'
                  pattern = sets.pvalue.pattern;
                case 'ecodam'
                  pattern = 'GATC';
              end
              randSeq = randseq(lenShortBp);
              [~, refDotsBpRes] = restrict( ...
                randSeq(:)', ...
                pattern(:)', ...
                round(length(pattern)/2));
              randShort = zeros(1, lenShortBp);
              refDotsBpRes = refDotsBpRes(2:end);
              randShort(refDotsBpRes) = bppx;
              randShort = convolve_bar(randShort, psfBp, lenShortBp);
              randShort = interp1(randShort, linspace(1, lenShortBp, lenShort));
          end
          % Perform comparison
          import DL.Hca.compute_random_scores
          [bestScoreVec(i), scoreMeanVec(i)] = compute_random_scores( ...
            randLong, ...
            randShort, ...
            stretchMax, ...
            sets.pvalue.barcodeType);
          % Update progress bar
          estTimeRem = round(toc/i*(sets.pvalue.numRnd-i)/6)/10;
          waitbar( ...
            i/sets.pvalue.numRnd, ...
            f, ...
            compose("Comparing barcode pair %.0f, estimated time remaining: %.1f minutes.", i+1, estTimeRem), ...
            'Position', [0, 0, 200, 50]);
        end
        close(f)
        data{end+1} = bestScoreVec;
        means{end+1} = scoreMeanVec;
        % Export best scores
        import CBT.Hca.Export.export_pval_struct;
        export_pval_struct(scoreFilePath, vals, data);
        export_pval_struct(meansFilePath, vals, means);
      end
    end
  end
end