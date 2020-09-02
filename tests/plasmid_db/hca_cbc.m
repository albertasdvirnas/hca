function [] = hca_plasmid_db_consensus()
    % run through the plasmid database to compute consensus for each
    % theory. Important to see what is the maximum PCC. how does it relate
    % to how many timeframes we take. Also how many barcodes we take.
    
    resultFold = 'resultData';
    setsF = 'sets_plasmid_db.txt';
    
    % link to plasmid db
    source = '/media/albyback/My Passport/DATA/plasmidDB/';
    
    % create tifstorunfile for each folder
    d = dir(source);
    dfolders = d([d(:).isdir]==1);

    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

    tifsName = cell(1,length(dfolders));
    
    for idx=1:length(dfolders)
        fileList = dir(fullfile(dfolders(idx).folder,dfolders(idx).name, 'Raw Kymographs','*.tif'));
        fileList = {fileList(:).name}';
        tifsName{idx} = fullfile(resultFold,strcat(['tifs_' dfolders(idx).name '.txt']));
        fid = fopen(tifsName{idx},'w');
        fprintf(fid,'%s\n',fileList{:});
        fclose(fid);
    end
    
    sets= ini2struct(setsF);
    sets.timeFramesNr = 0;

     % now for this, we run
    consensusStructs = cell(1,length(tifsName));
    kymoStructs = cell(1,length(tifsName));
    barcodeGen = cell(1,length(tifsName));

     for idx =1:length(tifsName)
        [consensusStructs{idx}, kymoStructs{idx},barcodeGen{idx}] = hca_no_theory( tifsName{idx}, sets );
     end
     
     idx = 1;
%      hcaSessionStruct
    sets.output = 'resultData';
    if ~exist(sets.output )
        mkdir(sets.output);
    end
    
    import CBT.Hca.Export.export_cbc_compatible_consensus;
    export_cbc_compatible_consensus(consensusStructs{idx}, barcodeGen{idx},kymoStructs{idx},sets);

    %% concentric plots for all: Note that we have function for plotting best concentric image in
%         import CBT.Hca.UI.Helper.plot_best_concentric_image;
%     plot_best_concentric_image(hAxis,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets);

    idx = 5;


          % loads figure window
    hFig = figure(...
        'Name', 'CB HCA tool', ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none' ...
    );
    hMenuParent = hFig;
    hPanel = uipanel('Parent', hFig);
    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);
    
    idx2 = [1:length(consensusStructs{idx}.consensusStruct)];
    import CBT.Hca.Export.plot_consensus_concentrically;
    plot_consensus_concentrically(consensusStructs{idx},barcodeGen{idx},idx2)
    
    
    
    import CBT.Hca.Core.gen_consensus_plot;
    gen_consensus_plot(consensusStructs,sets);
    
    
    % want to make a CBC plot like in CBC, but for Hca barcodes..
  

    % want to make this into full supported CBC 
                

%     % now for this, we run
%     consensusStructs = cell(1,length(tifsName));
%     kymoStructs = cell(1,length(tifsName));
%     barcodeGen = cell(1,length(tifsName));
% 
% 	 % load settings
%     sets = ini2struct( setsF );
%     idy = 1:5:100
%     maxCC = zeros(length(idy),length(tifsName));
%     for j=1:length(idy)
%         j
%         sets.timeFramesNr = idy(j);
% 
%         for idx=1:length(tifsName)
% 
%             [consensusStructs{idx}, kymoStructs{idx},barcodeGen{idx}] = hca_no_theory( tifsName{idx}, sets );
%         end
%         maxCC(j,:) = cellfun(@(x) x.treeStruct.maxCorCoef(1),consensusStructs);
%     end
        

end

