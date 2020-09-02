function [ sets ] = get_user_settings( sets )
    % get_user_settings
    %
    % This function asks for all the required settings for the generation
    % of the results at this point.
    %
	%     Args:
    %         sets (struct): Input settings to the method
    % 
    %     Returns:
    %         sets: Return structure

    if ~sets.kymosets.askforkymos 
        try 
            fid = fopen(sets.kymosets.kymoFile); 
            fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
            for i=1:length(fastaNames{1})
                [FILEPATH,NAME,EXT] = fileparts(fastaNames{1}{i});

                sets.kymosets.filenames{i} = strcat(NAME,EXT);
                sets.kymosets.kymofilefold{i} = FILEPATH;
            end
        catch
            sets.kymosets.askforkymos   = 1;
        end
%         sets.kymosets.kymoFile = 'kymos_2020-03-13_14_06_31_.txt';
    end
    
    if sets.kymosets.askforkymos 
        % loads figure window
        import Fancy.UI.Templates.create_figure_window;
        [hMenuParent, tsHCA] = create_figure_window('CB HCA tool','HCA');

        import Fancy.UI.Templates.create_import_tab;
        cache = create_import_tab(hMenuParent,tsHCA,'kymo');
        uiwait(gcf);  
        
        dd = cache('selectedItems');
        sets.kymosets.filenames = dd(1:end/2);
        sets.kymosets.kymofilefold = dd((end/2+1):end);
        delete(hMenuParent);
    end
    
    % add kymograph settings
	if sets.kymosets.askforsets
        prompt = {'Number of time frames (all time frames by default)','Alignment method (1 is nralign, 2 is ssdalign)', 'Filter the barcodes', 'Add consensus','Non-default bitmask settings', 'Non-defaut edge detection settings', 'Skip edge detection', 'Generate random fragments cut-outs','Generate independent subfragments for barcodes','Comparison Method'};
        title = 'Kymograph settings';
        dims = [1 35];
        definput = {'0','1','0','0','0','0','0','0','0','mass_pcc'};
        answer = inputdlg(prompt,title,dims,definput);
        
        sets.timeFramesNr = str2double(answer{1});
        sets.alignMethod = str2double(answer{2});
        sets.filterSettings.filter = str2double(answer{3});
        sets.genConsensus = str2double(answer{4});
        sets.bitmaskSettings = str2double(answer{5});
        sets.edgeSettings = str2double(answer{6});
        sets.skipEdgeDetection = str2double(answer{7}); % skips edge detection and assumes the first
        % non-zero and last nonzero to be edges of molecule
        sets.random.generate = str2double(answer{8});
        sets.subfragment.generate = str2double(answer{9});
        sets.comparisonMethod = answer{10};
    
        if sets.edgeSettings == 0
            sets.skipDoubleTanhAdjustment = 1;
        end
        
        % add filter settings if filtering is enabled
        if sets.filterSettings.filter
            % here we add settings
            prompt = {'Filter method','Filter size (px)'};
            title = 'Filtering settings';
            dims = [1 35];
            definput = {'0','2.3'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.filterSettings.filterMethod = str2double(answer{1}); % filter before or after stretching
            sets.filterSettings.filterSize = str2double(answer{2});      
        end

        if sets.genConsensus
           % here we add settings for this 
            prompt = {'Averaging method','Prompt for threshold', 'Consensus threshold'};
            title = 'Consensus settings';
            dims = [1 35];
            definput = {'bgmean','1','0.75'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.consensus.barcodeNormalization = answer{1};
            sets.consensus.promptForBarcodeClusterLimit =  str2double(answer{2}); 

            sets.consensus.threshold =  str2double(answer{3}); 
        end

        if sets.bitmaskSettings
            % here we add settings
            prompt = {'Camera (nm/px)','Point spread function (nm)', 'Delta cut'};
            title = 'Bitmasking settings';
            dims = [1 35];
            definput = {'130','300','3'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.bitmasking.prestretchPixelWidth_nm = str2double(answer{1}); % camera nm/px
            sets.bitmasking.psfSigmaWidth_nm = str2double(answer{2}); % psf in nanometers
            sets.bitmasking.deltaCut = str2double(answer{3}); % how many delta's should we take

            % number of untrusted pixels
            sets.bitmasking.untrustedPx = sets.bitmasking.deltaCut*sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
        end
    
        % edge settings. Add a new possibility, of just fitting the sigmoid
        % with new method introduced at v0.4
        if sets.edgeSettings
            answer = questdlg('Which edge detection method you want to choose', ...
            'Edge detection settings', ...
            'Otsu','Double tanh','Error function','Otsu');
            % Handle response
            switch answer
                case 'Otsu'
                    sets.edgeDetectionSettings.method = 'Otsu';
                case 'Double tanh'
                    sets.edgeDetectionSettings.method = 'Double tanh';
                case 'Error function'
                    sets.edgeDetectionSettings.method = 'Error Function';
            end

            
            
%             % add edge detection choices here
%             prompt = {'Skip double tan (default is 0, use double tan)'};
%             title = 'Edge detection settings';
%             dims = [1 35];
%             definput = {'0'};
%             answer = inputdlg(prompt,title,dims,definput);

%             sets.skipDoubleTanhAdjustment = str2double(answer{1});

            % would like to adjust these based on the experimental conditions!
            import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
            sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
        end

        if sets.random.generate
            % here we add settings
            prompt = {'Number of cut-outs','Length of a cut-out'};
            title = 'Cut-out settings';
            dims = [1 35];
            definput = {'10','50'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.random.noOfCutouts = str2double(answer{1}); % number of random cutouts from the input set
            sets.random.cutoutSize = str2double(answer{2});    
        end

        if  sets.subfragment.generate 
            prompt = {'Number of subfragments'};
            title = 'Number of subfragments settings';
            dims = [1 35];
            definput = {'2'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.subfragment.numberFragments = str2double(answer{1}); % number of random cutouts from the input set
        end
        
        if sets.output.askforoutputdir
            sets.output.matDirpath = strcat([uigetdir(pwd,'Choose a folder where you want to save the output') '/']);
            mkdir(sets.output.matDirpath);
        end
    end
    
end

