function [consensusStructs, kymoStructs,barcodeGen] = hca_no_theory( tifs, sets )
  % hca_no_theory
   
    % settings - name of the settings file
    % tif - name of the tif file
    
    % no output displayed apart from messages to the terminal, the 
    % output saved in the output folder
    
    % In case input files are not provided, set them automatically
    % ???
    if nargin < 1
        tifs = 'tifstorun.txt';
    end
    

% 
%     import CBT.Hca.Settings.set_fast_sets;
%     sets = set_fast_sets();

    
%     % load settings to a structure
%     try fid=fopen(settings); C2 = textscan(fid, '%s %s'); fclose(fid);
%         sets = cell2struct(C2{2},C2{1},1);
%     catch
%         error('No valid settings file provided');
%     end
%     
%     
%     % convert names to numerical values, apart from the values which should
%     % be strings
%     fn = fieldnames(sets);
%     for k=1:numel(fn)
%         if( ~isequal(fn{k},'literature'))
%             sets.(fn{k}) = str2num(sets.(fn{k}));            
%         end
% %     end
%    try 
%         fid = fopen(theoryfiles); 
%         theories = textscan(fid,'%s','delimiter','\n'); fclose(fid);
%         for i=1:length(theories{1})
%             [sets.theoryFileFold{i}, name, ext] = fileparts(theories{1}{i});
%             sets.theoryFile{i} = strcat([name ext]);
%         end
%    catch
%         error('No valid theories provided, please check the provided file');
%    end
%    
%    
 
   try 
        fid = fopen(tifs); 
        tifNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
   catch
        error('No valid tifs provided, please check the provided file');
   end
    
    kymoStructs = cell(1,length(tifNames{1}));
    
    % loop over movie file folder
    for idx = 1:length(tifNames{1})     
        try
            kymoStructs{idx}.unalignedKymo = imread(tifNames{1}{idx});
        catch
            fid = fopen(tifNames{1}{idx}); 
            file = textscan(fid,'%f '); fclose(fid);
            kymoStructs{idx}.unalignedKymo = file{1}';
        end
        kymoStructs{idx}.name = tifNames{1}{idx};
    end


%       % add kymographs
%     import CBT.Hca.Import.add_kymographs_fun;
%     [kymoStructs] = add_kymographs_fun(sets);

    %  put the kymographs into the structure
    import CBT.Hca.Core.edit_kymographs_fun;
    kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

    % align kymos
    import CBT.Hca.Core.align_kymos;
    [kymoStructs] = align_kymos(sets,kymoStructs);
% 
%     for i=1:length(kymoStructs)
%         kymoStructs{i}.alignedKymo = kymoStructs{i}.alignedKymo';
%     end
    % generate barcodes
    import CBT.Hca.Core.gen_barcodes;
    barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);

    
    % save kymos
    if sets.kymosets.savekymos
        currDate = strrep(datestr(datetime), ':', '_');
        mkdir('resultData',currDate); % makefolder with date
        for i =1:length(barcodeGen)

            % save current theory in txt file
            C = strsplit(barcodeGen{i}.name,'/');
            matFilename2 = strcat(['timeseries_' erase(C{end},'_kymograph.tif') '.txt']);
            matFilepath = fullfile('resultData', currDate,matFilename2);
            fd = fopen(matFilepath,'w');
            fprintf(fd, '%5.16f ', barcodeGen{i}.rawBarcode(barcodeGen{i}.rawBitmask));
            fclose(fd);

        end
    end
        
    
    
%     % generate consensus
%     import CBT.Hca.Core.gen_consensus;
%     consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);
    %

    % select consensus
%     import CBT.Hca.UI.Helper.select_consensus;
%     [consensusStruct,sets] = select_consensus(consensusStructs,sets);

% 
%     % generate random cut-outs
%     if sets.random.generate
%         % % Here generate cutouts using M_files_HCA_exp_cutouts:
%         import CBT.Hca.Core.Random.cutout_barcodes;
%         barcodeGenRandom = cutout_barcodes(barcodeGen,sets);
%         barcodeGenC = barcodeGenRandom;
%     else
%         barcodeGenC = barcodeGen;
%     end     
%     
% %     fname = strcat([sets.theoryFileFold{idx} 'rawBarcode.txt']);
% %     fileID = fopen(fname,'w');
% %     fprintf(fileID,'%2.5f ',barcodeGenC{1}.rawBarcode);
% %     fclose(fileID);
%          
% %     sets.theoryFile=[];
% %     sets.theoryFileFold = [];
% 
%     
% 	% now load theory
%     import CBT.Hca.UI.Helper.load_theory;
%     theoryStruct = load_theory(sets);
%     
% 	import CBT.Hca.Core.Analysis.convert_nm_ratio;
%     theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct );
%  
%     sets.comparison.useDTW = 0;
   % compare theory to experiment
%     import CBT.Hca.Core.Comparison.compare_theory_to_exp;
%     comparisonStructAll = compare_theory_to_exp(barcodeGenC, theoryStruct, sets, consensusStruct);
% 
%     
%     % bugcheck: if only one theory
%     import CBT.Hca.UI.combine_chromosome_results;
%     comparisonStruct = combine_chromosome_results(theoryStruct,comparisonStructAll);
% 
%     sets.output.matDirpath = 'resultData/';
%     
%     import CBT.Hca.Core.additional_computations
%     additional_computations(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct,comparisonStructAll,sets );
        
           %
  %  sets.displayResults = 1;
%     import CBT.Hca.UI.get_display_results;
%     get_display_results(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, sets);

%     
%     barNr = 1;
%     theoryNr = 8;
%     import CBT.Hca.Export.export_exp_vs_theory;
%     export_exp_vs_theory(comparisonStructAll,theoryStruct,barcodeGenC,sets,barNr,theoryNr);


%         
%         
%     theoryGen.sets = sets.theoryGen;
%     
%  	timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% 	matFilename = strcat(['theoryStruct_' sprintf('%s_%s', timestamp) 'session.mat']);
%     matFilepath = fullfile('resultData', matFilename);
% 
%     save(matFilepath, 'theoryGen');
%     
%     fprintf('Saved theory struct data ''%s'' to ''%s''\n', matFilename, matFilepath);


end

