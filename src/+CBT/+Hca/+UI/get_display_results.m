function [] = get_display_results(barcodeGen, consensusStruct,comparisonStruct,theoryStruct, sets)
    % get_display_results 
    % Display results from the comparison of experiments vs theory
    %     Args:
    %         barcodeGen: Input 
    %         consensusStruct: 
    %         comparisonStruct: 
    %         theoryStruct:
    %         sets:

    % if it was chosen to display results
    if sets.displayResults ==0
        fig1 = figure('Visible', 'off');
    else
        fig1 = figure;
    end

    % choose markers for everything
    markers = ['o';'s';'x';'+';'d';'v'];

    % compute cummulative sum of lengths of barcodes
    lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));

    % how many barcodes were compared to
    numBar = length(comparisonStruct);

    % If consensus was generated, then actual number of barcodes is one less
    if sets.genConsensus == 1
        numBar = numBar-length(consensusStruct);
    end

    % plot max corr coefs
    subplot(2,2,1);hold on;
    import CBT.Hca.UI.Helper.plot_max_coef;
    [fig1,maxcoef] = plot_max_coef(fig1,comparisonStruct, numBar, sets, markers);

    % plot best positions
    subplot(2,2,2);hold on;
    import CBT.Hca.UI.Helper.plot_best_pos;
    plot_best_pos(fig1,comparisonStruct, numBar, sets, markers,lengthBorders);

    subplot(2,2,3), hold on
%     if isequal(sets.comparisonMethod,'dtw')
%         
%     else
    import CBT.Hca.UI.Helper.plot_best_bar;
    plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef);
%     end

    

            
%     subplot(2,2,4), hold on
% %     option: alternatively plot concentirc plot of the two here, based on
%     % user input
%     import CBT.Hca.UI.Helper.plot_best_image;
%     plot_best_image(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef);
% %     fig1=figure;
        hAxis = subplot(2,2,4); hold on

    import CBT.Hca.UI.Helper.plot_best_concentric_image;
    plot_best_concentric_image(hAxis,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets);

    % option: plot all possible matches, i.e. for all 
    disp( strcat(['Number of timeframes for the unfiltered barcodes were = ' num2str(sets.timeFramesNr)]));
    
    
    
    %% additional things/options, added in HCA 4.1

    % take a simple example for concentric plot
%     idx = 87;
  
        
    % todo: make more user friendly..
    try
        mkdir(sets.output.matDirpath,sets.timestamp);
    end
    
    try
        saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'result_plot.eps'),'epsc');
    end
       
    try
        if sets.plotallmatches == 1
             mkdir(fullfile(sets.output.matDirpath,sets.timestamp),'Plots');
            for i=1:size(maxcoef,1)
                max2 = nan(size(maxcoef));      max2(i,1) = maxcoef(i,1);
                fig1 = figure('Visible', 'off');
%                 fig1 = figure;
                plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, max2);
                saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([num2str(i) '_plot.jpg'])));

%                 saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([sets.timestamp '_' num2str(i) '_plot.eps'])),'epsc');

            end
        end
    end
    %    assignin('base','hcaSessionStruct',hcaSessionStruct)
    
  %  cache('hcaSessionStruct') = hcaSessionStruct ;
%     end
end