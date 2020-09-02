function [hcaStruct] = run_parallel_from_barcodes(setsFile, theories, barcodeGenC, numFrames, nmBp,comparisonMethod)
    % HCA_Gui
    % Used for comparing fagments of human chromosome to chromosomes using CB (competitive binding) theory   
    %
    %     Args:
    %         sets (struct): Input settings to the method
    %         hcaStruct (struct): Input structure, if non-empty, load
    %         result structure instead of computing everything
    % 
    %     Returns:
    %         hcaStruct: Return structure
    % 
    %     Example:
    %         This is an example: run [hcaStruct] = HCA_Gui(sets) 
    %
    
    % TODO: return warnings to the places where it is important to know
    % what is being done and to check for consistency
    
    % timestamp for the results
    
    import CBT.Hca.Import.import_hca_settings;

    if nargin < 1 % if settings were not provided
        % import settings
        setsFile = 'hca_parallel_settings.txt';
    end

    [sets] = import_hca_settings(setsFile);
    
        
    if nargin >= 2 % if settings were not provided
    % import settings
        sets.theories = theories;
    end
    
%             
%     if nargin >= 3 % if settings were not provided
%     % import settings
%         listing = dir(fullfile(kymoFile,'*.tif'));
% 
%         fd =fopen('tifs.txt','w');
%         for i=1:length(listing)
%             fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
%         end
%         fclose(fd);
% 
%         sets.kymosets.kymoFile = 'tifs.txt';
%     end
%     
        if nargin >= 4 % if settings were not provided
    % 
        sets.timeFramesNr = numFrames;
        end
        
        if nargin >= 5
             sets.theory.nmbp = nmBp;
        end
        
        if nargin >= 6
             sets.comparisonMethod = comparisonMethod;
        end
    
        consensusStruct = [];
        %% now user theories. They could already be in txt files (if generated 
        % (with HCA 4.0.0), but we should keep support for older theory files too 
        sets.theoryFile=[];
        sets.theoryFileFold = [];

        % get user theory
        import CBT.Hca.Settings.get_user_theory;
        [theoryStruct, sets] = get_user_theory(sets);

        % maybe want to redo these!
        
        %         % compare theory to experiment
        import CBT.Hca.Core.Comparison.compare_distance;
        [rezMax,bestBarStretch,bestLength] = compare_distance(barcodeGenC, theoryStruct, sets, consensusStruct);

        comparisonStructAll = rezMax;
        for i=1:length(comparisonStructAll)
            for j=1:length(bestBarStretch{i})
                comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
                comparisonStructAll{i}{j}.length = bestLength{i}(j);
            end
        end

%         % compare theory to experiment
%         import CBT.Hca.Core.Comparison.compare_theory_to_exp;
%         comparisonStructAll = compare_theory_to_exp(barcodeGenC, theoryStruct, sets, consensusStruct);
        import CBT.Hca.Core.Comparison.combine_theory_results;
        [comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

% combine_theory_results
%         % bugcheck: if only one theory
%         import CBT.Hca.UI.combine_chromosome_results;
%         comparisonStruct = combine_chromosome_results(theoryStruct,rezMax,bestBarStretch,bestLength);


   % assign all to base
%         assignin('base','barcodeGenC', barcodeGenC)
%         assignin('base','theoryStruct', theoryStruct)
%         assignin('base','comparisonStructAll', comparisonStructAll)
%         assignin('base','sets', sets)
%         assignin('base','comparisonStruct', comparisonStruct)
%         assignin('base','consensusStruct', consensusStruct)

        %
%        sets.displayResults = 1;
%         import CBT.Hca.UI.get_display_results;
%         get_display_results(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, sets);
% 
% 
%         import CBT.Hca.Core.additional_computations
%         additional_computations(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct,comparisonStructAll,sets );
%         
%         % extra test: check if barcodes are placed correctly where they
%         % should be placed (or approximately around the correct place)
%         % also: if they are stretched as they should be stretched
% 
% %         allowedErrorPos,allowedErrorStr
%         import CBT.Hca.Core.Comparison.find_correct_matchings;
%         [tpr,tp,pStr] = find_correct_matchings(theoryStruct,barcodeGenC,comparisonStruct);
%         import CBT.Hca.UI.compute_true_positives;


        hcaStruct.barcodeGenC = barcodeGenC;
        hcaStruct.theoryStruct = theoryStruct;
        hcaStruct.comparisonStructAll = comparisonStructAll;
        hcaStruct.sets = sets;
        hcaStruct.comparisonStruct = comparisonStruct;
%         hcaStruct.consensusStruct = consensusStruct;
%         hcaStruct.consensusStructs = consensusStructs;
%         hcaStruct.kymoStructs = kymoStructs;

      %  sets.output.matDirpath = '/home/albyback/rawData/dnaData/humanData/output/';
  
%     else
%         % run load results function or something similar..
%     end
    
end