function [theoryGen] = hca_theory_script( fastas,settings )
  % hca_theory_script
  % This script generates the theory of selected fasta files
  % mimicing the hca theory gui, and runnable from the terminal
    
    % settings - name of the settings file
    % tif - name of the tif file
    
    % no output displayed apart from messages to the terminal, the 
    % output saved in the output folder
    
    % In case input files are not provided, set them automatically
    if nargin < 2
        settings = 'sets.txt';
        fastas = 'fastatorun.txt';
    end
    
    sets = ini2struct( settings );
    
    import CBT.Hca.Core.Theory.cb_model;
    sets.model = cb_model();



    % compute free concentrations
    import CBT.Hca.Core.Theory.compute_free_conc;
    sets = compute_free_conc(sets);

    theoryGen = struct();

   try 
        fid = fopen(fastas); 
        fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
   catch
        error('No valid fasta provided, please check the provided file');
   end
    
    import CBT.Hca.Core.Theory.create_memory_struct;
    import CBT.Hca.Core.Theory.compute_hca_theory_barcode;
    import CBT.Hca.Core.Theory.compute_theory_barcode;
    
    mkdir('resultData');
    % loop over movie file folder
    for idx = 1:length(fastaNames{1})
        %addpath(genpath(fastaNames{idx}))
        
        [chr1,header] = create_memory_struct(fastaNames{1}{idx});
        %
%         chr1 = memmapfile('chr1.mm', 'format', 'uint8');        

        disp(strcat(['loaded theory sequence ' fastaNames{1}{idx}] ));
        
        % add possibility for alternative theory
        % split the theory into shorter fragments here already
          theorySeq = compute_theory_barcode(chr1,sets);
%         theorySeq = compute_hca_theory_barcode(chr1,sets);
%         import CBT.Hca.Core.Theory.compute_hca_theory_fast;
%         theorySeq = compute_hca_theory_fast(chr1,sets);
%         
        theoryGen.theoryBarcodes{idx} = theorySeq;
        theoryGen.theoryNames{idx} = header(2:end);
        theoryGen.theoryIdx{idx} = idx;
        theoryGen.bpNm{idx} = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;
        
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(sets.theoryGen.meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile('resultData', matFilename2);
        fd = fopen(matFilepath,'w');
        fprintf(fd, '%5.3f ', theorySeq);
        fclose(fd);
        
        delete(chr1.Filename)
    end

    theoryGen.sets = sets.theoryGen;
    
 	timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
	matFilename = strcat(['theoryStruct_' sprintf('%s_%s', timestamp) 'session.mat']);
    matFilepath = fullfile('resultData', matFilename);

    save(matFilepath, 'theoryGen');
    
    fprintf('Saved theory struct data ''%s'' to ''%s''\n', matFilename, matFilepath);


end

