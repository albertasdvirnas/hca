function [data, means] = run_random_score_gen(sets)
% run_random_score_gen
import CBT.Hca.Core.Pvalue.convolve_bar;

[outputFolder, ~] = fileparts(sets.pvalue.file);
if not(isfolder(outputFolder))
  mkdir(outputFolder);
end
if not(isfolder(sets.pvalue.fastaFolder))
  mkdir(sets.pvalue.fastaFolder);
end

% Variables
bppx = sets.pvalue.pixelWidth_nm/sets.pvalue.nmbp;
recSite = sets.pvalue.pattern;
lenLong_all = sets.pvalue.lenLong;
lenShort_all = sets.pvalue.lenShort;
psfPx_all = sets.pvalue.psfSigmaWidth_nm/sets.pvalue.pixelWidth_nm;
stretchMax_all = sets.pvalue.stretchMax;
switch sets.pvalue.variableType
  case 'lenLong'
    lenLong_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    vals = lenLong_all;
  case 'lenShort'
    lenShort_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    vals = lenShort_all;
  case 'psf'
    psfPx_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    vals = psfPx_all;
  case 'stretchMax'
    stretchMax_all = sets.pvalue.varMin:sets.pvalue.varStep:sets.pvalue.varMax;
    vals = stretchMax_all;
  otherwise
    error('Bad variable type selected')
end

% Output
import CBT.Hca.Import.load_pval_struct;
[~, data] = load_pval_struct(sets.pvalue.file);
[~, means] = load_pval_struct(sets.pvalue.meansFile);
alreadyInData = length(data);

% For each psf
for psfPx=psfPx_all
  psfBp = round(psfPx*bppx);
  % For each long length
  for lenLong=lenLong_all
    % Generate long barcode
    switch sets.pvalue.barcodeType
      case 'cb'
        randLong = normrnd(0, 1, 1, lenLong);
        randLong = convolve_bar(randLong, psfPx, lenLong);
      case 'dots'
        lenLongBp = round(lenLong*bppx);
        randFastaName = compose("random_fasta_length_%.0f.fasta", lenLongBp);
        randFastaPath = fullfile(sets.pvalue.fastaFolder, randFastaName);
        if not(isfile(randFastaPath))
          import DL.Hca.generate_random_fasta
          generate_random_fasta(randFastaPath, lenLongBp)
        end
        import DL.Hca.compute_theory_dotbar
        randLong = compute_theory_dotbar( ...
          randFastaPath, ...
          recSite, ...
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
        timespent = 0;
        % For each random sample
        for i=1:sets.pvalue.numRnd
          tic
          % Generate short barcode
          switch sets.pvalue.barcodeType
            case 'cb'
              randShort = normrnd(0, 1, 1, lenShort);
              randShort = convolve_bar(randShort, psfPx, lenShort);
            case 'dots'
              lenShortBp = round(lenShort*bppx);
              randSeq = randseq(lenShortBp);
              [~, refDotsBpRes] = restrict( ...
                randSeq(:)', ...
                recSite(:)', ...
                round(length(recSite)/2));
              randShort = zeros(1, lenShortBp);
              randShort(refDotsBpRes(2:end)) = bppx;
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
          timespent = timespent + toc;
          estTimeRem = round(timespent/i*(sets.pvalue.numRnd-i)/6)/10;
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
        export_pval_struct(sets.pvalue.file, vals, data);
        export_pval_struct(sets.pvalue.meansFile, vals, means);
      end
    end
  end
end
