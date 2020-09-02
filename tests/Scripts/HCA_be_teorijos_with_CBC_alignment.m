function HCA_be_teorijos_with_CBC_alignment(sets, hcaStruct)
    % HCA_be_teorijos
    % Used for comparing fagments of human chromosome to chromosomes using CB (competitive binding) theory 
    %
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
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

    if nargin < 1 % if settings were not provided
        % import settings
        import CBT.Hca.Import.import_hca_settings;
        [sets] = import_hca_settings('hca_settings.txt');
    end
    
    if nargin < 2 
        % load all user selected settings
        import CBT.Hca.Settings.get_user_settings;
        sets = get_user_settings(sets);
        
        % add kymographs
        import CBT.Hca.Import.add_kymographs_fun;
        [kymoStructs] = add_kymographs_fun(sets);

        %  put the kymographs into the structure
        import CBT.Hca.Core.edit_kymographs_fun;
        kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

        import CBT.Hca.Core.align_kymos_consensus;
        [ kymoStructs,barcodeGen, consensusStructs ] = align_kymos_consensus( sets, kymoStructs );
% 
%         import CBT.Hca.Core.align_kymos_consensus;
%         barcodeGen =  align_kymos_consensus(kymoStructs, sets);
%         
        
        % if we want to plot
        import CBT.Hca.Export.plot_consensus_concentrically;
        plot_consensus_concentrically(consensusStructs,barcodeGen)


        %
%         % align kymos
%         import CBT.Hca.Core.align_kymos;
%         [kymoStructs] = align_kymos(sets,kymoStructs);
%         
%         % generate barcodes
%         import CBT.Hca.Core.gen_barcodes;
%         barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);
%  
        % generate consensus
%         import CBT.Hca.Core.gen_consensus;
%         consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);
%         %

%         % select consensus
%         import CBT.Hca.UI.Helper.select_consensus
%         [consensusStruct,sets] = select_consensus(consensusStructs,sets);
        
        % TODO : add consensus plot
%         % recreate consensus result from cbc_gui? 
%         import CBT.Hca.Core.gen_consensus_plot;
%         gen_consensus_plot(consensusStructs,sets);
        
%         % generate random cut-outs
%         if sets.random.generate
%             % % Here generate cutouts using M_files_HCA_exp_cutouts:
%             import CBT.Hca.Core.Random.cutout_barcodes;
%             barcodeGenRandom = cutout_barcodes(barcodeGen,sets);
%             barcodeGenC = barcodeGenRandom;
%         else
%             barcodeGenC = barcodeGen;
%         end
% 
%     
%         %% now user theories. They could already be in txt files (if generated 
%         % (with HCA 4.0.0), but we should keep support for older theory files too 
%         sets.theoryFile=[];
%         sets.theoryFileFold = [];
% 
%         % get user theory
%         import CBT.Hca.Settings.get_user_theory;
%         [theoryStruct, sets] = get_user_theory(sets);
% 
%         % compare theory to experiment
%         import CBT.Hca.Core.Comparison.compare_theory_to_exp;
%         comparisonStructAll = compare_theory_to_exp(barcodeGenC, theoryStruct, sets, consensusStruct);
% 
% 
%         % bugcheck: if only one theory
%         import CBT.Hca.UI.combine_chromosome_results;
%         comparisonStruct = combine_chromosome_results(theoryStruct,comparisonStructAll);
% 
% 
%    % assign all to base
        assignin('base','barcodeGen', barcodeGen)
        assignin('base','kymoStructs', kymoStructs)
        assignin('base','consensusStructs', consensusStructs)

%         assignin('base','theoryStruct', theoryStruct)
%         assignin('base','comparisonStructAll', comparisonStructAll)
%         assignin('base','sets', sets)
%         assignin('base','comparisonStruct', comparisonStruct)
%         assignin('base','consensusStruct', consensusStruct)
% 
%         %
%       %  sets.displayResults = 1;
%         import CBT.Hca.UI.get_display_results;
%         get_display_results(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, sets);
% 
% 
%         import CBT.Hca.Core.additional_computations
%         additional_computations(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct,comparisonStructAll,sets );
%         
% 
%       
%       %  sets.output.matDirpath = '/home/albyback/rawData/dnaData/humanData/output/';
%   
    else
        % run load results function or something similar..
    end
    
end