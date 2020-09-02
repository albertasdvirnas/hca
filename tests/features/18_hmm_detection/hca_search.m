classdef hca_search
    %HCA_SEARCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sets;
        kymoStructs;
        kymoStructsAligned;
        barcodeGen;
        theoryStruct;
        randKymos;
        names;
        posCuts;lengths;
        theoriesTxt;
        ver;
        tifsTxt;
    end
    
    methods
        function obj = hca_search(settingsTxt,theoriesTxt,tifsTxt,outputDir)
           % settings
            % this could also be loaded as a class
            import CBT.Hca.Import.import_hca_settings;
            [obj.sets] = import_hca_settings(settingsTxt);
            
           
            obj.sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
            obj.sets.kymosets.kymoFile = tifsTxt;


            import CBT.Hca.Settings.get_user_settings;
            obj.sets = get_user_settings(obj.sets);
            obj.sets.output.matDirpath = outputDir;
            
            % for now keep all kymos
            obj.sets.whichtokeep = length(obj.sets.kymosets.filenames);
            
            % make directory if it doesn't exist
            try mkdir(obj.sets.output.matDirpath);catch; end

            % add kymographs
            import CBT.Hca.Import.add_kymographs_fun;
            [obj.kymoStructs] = add_kymographs_fun(obj.sets);

            %  put the kymographs into the structure
            import CBT.Hca.Core.edit_kymographs_fun;
            obj.kymoStructs = edit_kymographs_fun(obj.kymoStructs,0);

            % to multibar - this is optional.
            import generate.kymo_to_multi_bar;
            obj.kymoStructs = kymo_to_multi_bar(obj.kymoStructs);
            % now convert this to many kymo of single frame

            % align kymos - should not take any time if it's just single frame
            import CBT.Hca.Core.align_kymos;
            [obj.kymoStructsAligned] = align_kymos(obj.sets,obj.kymoStructs);

               % generate barcodes
            import CBT.Hca.Core.gen_barcodes;
            obj.barcodeGen =  CBT.Hca.Core.gen_barcodes(obj.kymoStructsAligned, obj.sets);
            %  
            % trueStart = -ver+50;
            % startFound = cellfun(@(x) x.lE,barcodeGen);

            % %% !!!
            % % this shows if we estimate the start pixel correctly (based on PSF)
            % figure,plot(startFound);hold on;plot(trueStart); legend({'Found','True'})
            % figure,plot(startFound-trueStart)

            obj.sets.theories = theoriesTxt;
            % this should be in settings file
            obj.sets.theory.nmbp = 0.3;
            % get user theory
            import CBT.Hca.Settings.get_user_theory;
            [obj.theoryStruct, obj.sets] = get_user_theory(obj.sets);


        end
        
    end
    
    methods (Static)
        
        function obj = generate_sim(settingsTxt,numRows,outFold)

            obj.sets = ini2struct(settingsTxt);
            % % chech if the sequences have to be reproducible
            rng(1,'twister');

            % could put this to hca_search too
            obj.sets.outFold = outFold;
            [obj.randKymos,obj.names,obj.posCuts,obj.lengths,obj.theoriesTxt,obj.ver,obj.tifsTxt] = generate_sim_theory_single(obj.sets.lengthN,obj.sets.kernelsigma,obj.sets.outFold,obj.sets.lambda,obj.sets.pccScore,obj.sets.strF,numRows);
            % 
            obj.sets.kymosets.kymoFile = obj.tifsTxt;
        end
    end
end

