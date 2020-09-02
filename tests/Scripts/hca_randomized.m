function [] = hca_randomized( theoryfiles)
  % hca_theory_script
  % This script generates the theory of selected fasta files
  % mimicing the hca theory gui, and runnable from the terminal
    
    % settings - name of the settings file
    % tif - name of the tif file
    
    % no output displayed apart from messages to the terminal, the 
    % output saved in the output folder

    import CBT.Hca.Settings.set_fast_sets;
    sets = set_fast_sets();

    
    % In case input files are not provided, set them automatically
    if nargin < 2
%         settings = 'sets.txt';
%         theory = 'tifstorun.txt';
        theoryfiles = 'theorytorandomizefiles.txt';
        sets.random.cutoutSize = 683;
        sets.random.noOfCutouts = 50;
    end
    

    try 
        fid = fopen(theoryfiles); 
        theories = textscan(fid,'%s','delimiter','\n'); fclose(fid);
        for i=1:length(theories{1})
            [sets.theoryFileFold{i}, name, ext] = fileparts(theories{1}{i});
            sets.theoryFile{i} = strcat([name ext]);
        end
    catch
        error('No valid theories provided, please check the provided file');
    end

    % load theory
    import CBT.Hca.UI.Helper.load_theory;
    theoryStruct = load_theory(sets);
    
    % nkow we loop through the theories to generate sets of randomized
    % barcodes
    
    rng('default');
    rng(1);
    
    for i=1:length(theoryStruct)
        
        % first load the barcode
        fileID = fopen(theoryStruct{i}.filename,'r');
        formatSpec = '%f';
        seq = fscanf(fileID,formatSpec);
        fclose(fileID);
        [sk y] = strtok(theoryStruct{i}.name,' ');
        theory_dir = fullfile('resultData', sk);
        mkdir(theory_dir);
        theory_mat = fullfile('resultData', 'mat');
        mkdir(theory_mat);

        % put this theory to barcodeGen
        barcodeGen = [];
        barcodeGen{1}.rawBarcode = seq';
        barcodeGen{1}.rawBitmask = ones(1,length(seq));

        import CBT.Hca.Core.Random.cutout_barcodes;
        barcodeGenRandom = cutout_barcodes(barcodeGen,sets);
        
        
        % generate random barcode of the same length:
        % later take this from settings.
        psf = 300/110;

        import CBT.Hca.Core.Pvalue.convolve_bar;
        import SignalRegistration.unmasked_pcc_corr;

        
        for k=1:sets.random.noOfCutouts
    
            % compute the long random barcode
            rand2 = normrnd(0, 1, 1, sets.random.cutoutSize);
            % convolve with a Gaussian
            rand2 = zscore(convolve_bar(rand2, psf, length(rand2)));

            bar = zscore(barcodeGenRandom{k}.rawBarcode);
            fun = @(alpha) (max(max(unmasked_pcc_corr(bar*alpha+(1-alpha)*rand2, bar, barcodeGenRandom{k}.rawBitmask)))-0.8).^2;
            alpha = lsqnonlin(fun, 0.8);

            noisified_barcode = alpha*bar + (1-alpha)*rand2;
            noisified_barcode(~barcodeGenRandom{k}.rawBitmask)  = nan;

    %         % Since we want only one position, could potentially speed this up.
    %         xcorrs = unmasked_pcc_corr(barcodeGenRandom{1}.rawBarcode, barcodeGenRandom{1}.rawBarcode, barcodeGenRandom{1}.rawBitmask);
    %         

            noisified_barcode = noisified_barcode - min(noisified_barcode(:));
            noisified_barcode = noisified_barcode/max(noisified_barcode(:));
            barcodeGenRandom{k}.noiseBarcode = noisified_barcode;
            barcodeGenRandom{k}.alpha = alpha;
            barcodeGenRandom{k}.noise = rand2;
            matFilepath = fullfile('resultData',sk, strcat(['noisified_theory_barcode_nr_' num2str(k) '_alpha_' num2str(alpha) '_' genvarname(sk) '.tif']));

            % now we save this noisified_barcode as a .tif
            imwrite(noisified_barcode, matFilepath);
        end
        matFilepath = fullfile(theory_mat, strcat(['noisified_theory_barcode_'  genvarname(sk) '.mat']));

        save(matFilepath, 'barcodeGenRandom');

    end
%     