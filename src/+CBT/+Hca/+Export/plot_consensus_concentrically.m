function [  ] = plot_consensus_concentrically(consensusStructs,barcodeGen, idx, ts )
    % plot_consensus_concentrically
    %
    %   Args:
    %       ts, consensusStructs,consensusStruct,barcodeGen
    %
    %   No returns
    
    % which indices of consensus to plot
    if nargin < 3
        % if nothing selected, plot all of them
        idx = 1:length(consensusStructs.consensusStruct);
    end
    
    import Fancy.UI.FancyTabs.TabbedScreen;

    if nargin < 4
%                 % loads figure window
%         hFig = figure( 'Units', 'normalized', ...
%             'OuterPosition', [0 0 1 1], 'NumberTitle', 'off', 'MenuBar', 'none' );
%         hPanel = uipanel('Parent', hFig);
%         ts = TabbedScreen(hPanel);
    end
    
% 
%     hTabConcentric = ts.create_tab('Concentric Plots');
%     hPanelConcentric = uipanel('Parent', hTabConcentric);
%     tsConcentric = TabbedScreen(hPanelConcentric);

    clK = cell(1,length(idx));
    for i=1:length(idx)
        curidx = idx(i);
%         num2str(cell2mat(consensusStructs.treeStruct.clusteredBar{curidx}));
        vec = cell2mat(consensusStructs.treeStruct.clusteredBar{curidx});
        clK{i} = strjoin(arrayfun(@(x) num2str(x),vec,'UniformOutput',false),',');
    end 
    
    numBarcodesInClusters = cellfun(@(clK) sum(clK == ','), clK) + 1;
    [~, sortOrder] = sort(numBarcodesInClusters, 'descend');
    numClusters = length(clK);

    
    f = figure('Name','Consensus', 'Units', 'normalized', ...
            'OuterPosition', [0 0 1 1], 'NumberTitle', 'off', 'MenuBar', 'none' );
        
    tabgp = uitabgroup(f,'Units', 'normalized','Position',[.0 .0 1 1]);
    tab1 = uitab(tabgp,'Title','Concentric plot');
    tab2 = uitab(tabgp,'Title','Linear plot');
    tab3 = uitab(tabgp,'Title','CC scores');

    tabgp2 = uitabgroup(tab1,'Units', 'normalized','Position',[.0 .0 1 1]);

  
    for clusterNum = 1:numClusters
%          clusterNum = 1
        clusterIdx = sortOrder(clusterNum);
        ck = clK{clusterIdx};
        hTabCluster = uitab(tabgp2,'Title',sprintf('%d', clusterIdx));

        alignedBarcodes = consensusStructs.treeStruct.barcodes{clusterIdx};
        %alignedBitmasks = isnan(consensusStructs.treeStruct.barcodes{end});
        clusterBarcodeKeys = find(consensusStructs.treeStruct.barMatrix{clusterIdx});
        clusterBarcodeNames = cellfun(@(x) x.name, barcodeGen,'UniformOutput',false);
        clusterBarcodeAliases = clusterBarcodeNames(clusterBarcodeKeys);
        
        sanitizedBarcodesAligned = num2cell(alignedBarcodes,2);
        clusterBarcodeNamesAligned = strcat({'Aligned Barcode '}, cellfun(@(x) num2str(x),num2cell(clusterBarcodeKeys),'UniformOutput',false), {' ('}, clusterBarcodeAliases, {')'});

%         tabTitleAligned = ck;
%         titleAligned = {'Aligned Barcodes in Consensus Cluster (Concentric Plot)',ck};

%         tabgp3 = uitabgroup(hTabCluster,'Units', 'normalized','Position',[.0 .0 1 1]);
%         hPanelAligned = uitab(tabgp3,'Title','Concentric plot','Position',[0, 0.1, 0.6, 0.8]);
%         hTable = uitab(tabgp3,'Title','Concentric plot','Position',[0.6, 0, 0.4, 1]);

%         tab2 = uitab(tabgp,'Title','Linear plot');

%             
%         hTabAligned = tsCluster.create_tab(tabTitleAligned);
%         tsCluster.select_tab(hTabAligned);
%         hPanelAligned = uipanel('Parent', hTabAligned);
%         hAxisAligned = axes(...
%             'Units', 'normal', ...
%             'Position', [0, 0.1, 0.6, 0.8], ...
%             'Parent', hPanelAligned ...
%             );

        uitable(hTabCluster,...
            'Data', clusterBarcodeNamesAligned',...'ColumnName'
            'ColumnName', {'Barcode Names (Numbered from center outwards)'},...
            'Units', 'normal', 'Position', [0.7, 0.3, 0.25, 0.6]);
       try
            sanitizedBarcodes = [sanitizedBarcodesAligned; {consensusStructs.consensusStruct{clusterIdx}.rawBarcode}];
            sanitizedBarcodes = cellfun(@(v) (v - nanmean(v(:)))./nanstd(v(:)), sanitizedBarcodes, 'UniformOutput', false);


            hAxis = axes('parent', hTabCluster,'Position', [0, 0.1, 0.6, 0.8]);
            title(sprintf('%s',ck) )
            import Barcoding.Visualizing.plot_circular_barcodes_concentrically;
            plot_circular_barcodes_concentrically(hAxis,sanitizedBarcodes);
       catch
       end
%         title(hAxisAligned, titleAligned);

    end
    
    
        tabgp3 = uitabgroup(tab3,'Units', 'normalized','Position',[.0 .0 1 1]);
        hTabCluster3 = uitab(tabgp3,'Title','PCC scores for adding new value to a cluster');
        hAxis3 = axes('parent', hTabCluster3,'Position', [0.1, 0.1, 0.6, 0.8]);
        plot(consensusStructs.treeStruct.maxCorCoef)
        xlabel('Cluster');
        ylabel('CC score');
%         title(sprintf('%s',ck) )
   
  
end

