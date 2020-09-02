classdef hmm_search
    %HMM_SEARCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sets
    end
    
    methods
        function obj = hmm_search(settings)
            obj.sets = ini2struct( settings );

        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
%     
%     % synthetic run should be similar to mbio_run or puuh_run, but in this case
% % we know the correct placements, correct rescaling factors, etc.
% 
% % first generate data
% 
% % first need some settings
% 
% sets.fold = 'simData/';
% sets.svType = 'Insertion';
% sets.numSamples = 1;
% sets.kernelsigma = 2.3;
% sets.pccScore = 0.9;
% sets.circ = 0;
% sets.rand.randmethod = 'imgaussfilt';
% mkdir(sets.fold);
% inLength = 50;
% sets.length1 = 500;
% import Rand.generate_one_type_sv;
% [randStruct] = generate_one_type_sv(sets,inLength);

%%
%  settings = 'sensitivity_settings.txt' ;
% 
% % what are natural values for parameters?
% 
% 
% numBar = 1;
% % generate data
% sets.svList = 1*ones(1,numBar);
% 
% % psf sigma
% sets.kernelsigma = 2.3;
% sets.circ = 0;
% % generate sample linear structural variations
% import Rand.generate_linear_sv;
% [bar1, bar2,matchTable]  = arrayfun(@(x) generate_linear_sv(len1, lenVar, x,sets),sets.svList,'UniformOutput',false);

% 
% sets.length = length(bar2{1});
% import Rand.add_noise_to_barcode;
% bar2  = cellfun(@(x) add_noise_to_barcode(x, sets),bar2,'UniformOutput',false);



% now add random stretch factor (based on normrnd) to bar1 on range from -20 to 20
% randStr = num2cell(1 + normrnd(0,10,1,numBar)/100);
% 
% import Comparison.interpolate_data;
% % interpolate data
% barS = cellfun(@(x,y) interpolate_data(x,length(x), length(x)*y,sets),bar1,randStr,'UniformOutput',false);

%%
% 
% %% test which is plotted in paper
% len1 = 500;
% lenVar = 50;
% 
% % sets.svList = [1 2 3 4];
% sets.svList = [1];
% 
% sets.kernelsigma = 2.3;
% 
% % generate sample linear structural variations
% import Rand.generate_linear_sv;
% [bar1, bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(len1, lenVar, x,sets),sets.svList,'UniformOutput',false);
% 
% import Rand.save_rand_barcodes;
% names = save_rand_barcodes(bar1,bar2,matchTable,sets.fold);
%  
% % plot
% import Plot.alignment;
% sets.fold = 'out';
% alignment(matchTable,bar1{1},bar2{1},sets,'bla')
% 
% % ( res_table, bar1D, bar2D,sets,name)
% pairs1 = [1]; % odd barcodes
% pairs2 = [1]; % even barcodes

sett = 'mbio.txt';

         
% % this is ok.
% 
% sets.length = length(bar2{1});
% import Rand.add_noise_to_barcode;
% bar2  = cellfun(@(x) add_noise_to_barcode(x, sets),bar2,'UniformOutput',false);
% 
% 
% 
% % now add random stretch factor (based on normrnd) to bar1 on range from -20 to 20
% randStr = num2cell(1 + normrnd(0,10,1,numBar)/100);
% 
% import Comparison.interpolate_data;
% % interpolate data
% barS = cellfun(@(x,y) interpolate_data(x,length(x), length(x)*y,sets),bar1,randStr,'UniformOutput',false);
% 
% % now the matchtable does not exactly match! we save these bars here



    iiN = 100;
    sets = ini2struct( sett );

    % now we have to decide which comparisons we want to do

    %%
    sets.fold = strcat(num2str(iiN),'dataSim/');
    sets.foldResult = strcat(num2str(iiN),'resultData/');

    mkdir( sets.fold);
    mkdir( sets.foldResult);

%     res = cell(1,length(pairs1));
    comparisonStruct = cell(1,length(pairs1));
%         iscircularq = 0;
%     iscirculard = 0;
% %     iscircularq = isempty(find(hcaStruct{pairs1(ii)}.consensusStruct.rawBitmask==0));
% %     iscirculard = isempty(find(hcaStruct{pairs2(ii)}.consensusStruct.rawBitmask==0));
% 
%     % set if this comparison is circular (query needs to be circular
    sets.circ = 0;


for ii = 1:length(pairs1);
    disp(strcat(['Analyzing ' names{pairs1(ii)} ]));
    disp(strcat([names{pairs2(ii)} ]));

    % load data
    query_txt = bar1{pairs1(ii)};
    data_txt = bar2{pairs2(ii)};
    
    
    % how to deal with nan's. I think if there are nan's, then one of the
    % barcodes has to be treated as linear. Data is usually treaded as linear
    % (at least for this project). But maybe data should be treated as
    % circular? (Since data is usually a plasmid, 

    % choose only non nan's for comparison
%     query_txt = query_txt(hcaStruct{pairs1(ii)}.consensusStruct.rawBitmask);
%     data_txt = data_txt(hcaStruct{pairs2(ii)}.consensusStruct.rawBitmask);

    % first try to align query to data ( or oposite) using HCA (or maybe MP?), to find the
    % best stretching factor

    import mp.mp_determine_best_stretch;
    [bestStretch,maxCoef] = mp_determine_best_stretch(query_txt,data_txt,sets);
    
    % more accurate would be to run a new one for each stretch factor, but
    % we do just for the best one here.
    sets.lenVar = length(data_txt);
    import Core.run_pval_matlab;
    [thresh] = run_pval_matlab(sets.r, round(length(query_txt)*bestStretch),sets);
 

     if max(maxCoef) > thresh
         % then we continue

        % set if this comparison is circular
        sets.circ = iscircularq;

        import  Core.run_hmm_bash_stretch;
        sets.analysis.stretchF = 0.04;
        [resultStruct] = run_hmm_bash_stretch(query_txt, data_txt', sets, bestStretch, iiN);

        sets.lenVar = round(length(query_txt)*resultStruct.bestSt);
        
        % merge table
        % 
        
%         resultStruct
        if sets.merge
            import functions.merge_matchtable_full;
            [mergedMatchTable,barsA,barsB] = merge_matchtable_full(resultStruct,sets);
            res.matchTable = mergedMatchTable;
        else
            import functions.create_full_table;
            [~,barsA,barsB] = create_full_table( resultStruct.matchTable ,resultStruct.bar1,resultStruct.bar2, 1 );
             res.matchTable =   resultStruct.matchTable ;
        end
        
        res.bar1 = resultStruct.bar1;
        res.bar2 = resultStruct.bar2;
        res.fragmentpcc = cellfun(@(x,u) zscore(x,1)*zscore(u',1)/length(x), barsA,barsB);

        if sets.shift
            % now shift one of these. This should be tested via tests before
            % generating final figures for the paper
           shift = res.matchTable(1,1)-1;
           res.bar1  = circshift(res.bar1, [0,-shift]);
           for k=1:size(res.matchTable,1)
                res.matchTable(k,1:2) = mod( res.matchTable(k,1:2)-1-shift,length(res.bar1))+1; 
           end
        end
       res.lengths = cellfun(@(x) length(x),barsA);
        
        % might perform oddly if merge match table is used.. anyway only
        % for
        % visual reference
%         import Plot.plot_sv_one;
%         plot_sv_one( {res},sets,strcat(num2str(ii), 'puuh_plot.eps'),1,'title');
 
        import Core.run_pval_matlab;
        [thresh,par1,par2] = arrayfun(@(x) run_pval_matlab(x, length(data_txt),sets),res.lengths);
       
        passthresh = res.fragmentpcc > thresh;
        
        % compute p-vals
        pvalFun = @(x,y,z) 1-(0.5+vpa(0.5)*(vpa(1,16)-vpa(betainc(x.^2,0.5,y/2-1,'upper'),16))).^z;
        pval = arrayfun(@(x,y,z) pvalFun(x,y,z),res.fragmentpcc,par1,par2);

        % save to output
        comparisonStruct{ii}.raw.fragmentpcc = res.fragmentpcc;
        comparisonStruct{ii}.raw.matchTable = res.matchTable;
        comparisonStruct{ii}.raw.lengths = res.lengths;
        import Comparison.interpolate_data;
        comparisonStruct{ii}.raw.bar1 =   res.bar1 ;
        comparisonStruct{ii}.raw.bar2 = data_txt;
        comparisonStruct{ii}.raw.subseqA = barsA;
        comparisonStruct{ii}.raw.subseqB = barsB;
        comparisonStruct{ii}.raw.pass = passthresh;
        comparisonStruct{ii}.raw.cthresh = thresh;
        comparisonStruct{ii}.raw.pval = pval;


        comparisonStruct{ii}.filt.fragmentpcc =  res.fragmentpcc(passthresh);
        comparisonStruct{ii}.filt.matchTable = res.matchTable(passthresh,:);
        comparisonStruct{ii}.filt.lengths =res.lengths(passthresh);
        comparisonStruct{ii}.filt.bar1 = comparisonStruct{ii}.raw.bar1 ;
        comparisonStruct{ii}.filt.bar2 = data_txt;
        comparisonStruct{ii}.name1 = names{pairs1(ii)};
        comparisonStruct{ii}.name2 = names{pairs2(ii)};
     else
         comparisonStruct{ii} = [];
     end

end

    % want to shade the region of bars which do not match. Plot all
    % comparisons circularly? Add P-val!
    sets.fold = '/home/albyback/git/sv/manuscript/figs/';

    for ii=1:length(pairs1)
        d1=strsplit(comparisonStruct{ii}.name1,'/');
        d2=strsplit(comparisonStruct{ii}.name2,'/');
        if sets.merge
            import Plot.plot_sv_final;
            plot_sv_final( {comparisonStruct{ii}.raw},sets,strcat(num2str(ii), 'merge_puuh_plot.eps'),1,strcat([d1{5} '  vs. ' d2{5} ]));
        else
        import Plot.plot_sv_final;
        plot_sv_final( {comparisonStruct{ii}.raw},sets,strcat(num2str(ii), 'puuh_plot.eps'),1,strcat([d1{5} '  vs. ' d2{5} ]));
        end
    end
    % when plotting, shade the barcodes which do not pass the threshold
    
end

