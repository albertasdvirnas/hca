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

    % now for this, we run
    consensusStructs = cell(1,length(tifsName));
    kymoStructs = cell(1,length(tifsName));
    barcodeGen = cell(1,length(tifsName));

	 % load settings
    sets = ini2struct( setsF );
    idy = 1:5:100
    maxCC = zeros(length(idy),length(tifsName));
    for j=1:length(idy)
        j
        sets.timeFramesNr = idy(j);

        for idx=1:length(tifsName)

            [consensusStructs{idx}, kymoStructs{idx},barcodeGen{idx}] = hca_no_theory( tifsName{idx}, sets );
        end
        maxCC(j,:) = cellfun(@(x) x.treeStruct.maxCorCoef(1),consensusStructs);
    end
        

end

