function clusterConsensusData = export_old_cbc_format(kymoStructs,consensusStructs,consensusStruct,barcodeGen,sets)

    if nargin < 6
        idx = length(consensusStructs.treeStruct.maxCorCoef);
    end
    % now generate the output which is equivalent to CBC output (so it is
    % loadable with CBC software completely

    hcaClusterConsensusData = [];
    hcaClusterConsensusData.barcode = consensusStruct.rawBarcode;
    hcaClusterConsensusData.bitmask = consensusStruct.rawBitmask;
    hcaClusterConsensusData.indexWeights = consensusStruct.indexWeights;
    hcaClusterConsensusData.stdErrOfTheMean = consensusStruct.stdErrOfTheMean;
    hcaClusterConsensusData.clusterKey = num2str(sort(consensusStruct.indices));
    hcaClusterConsensusData.clusterKey = regexprep(hcaClusterConsensusData.clusterKey,' +',',');
    % strrep(,'   ', ',')
    % hcaClusterConsensusData.clusterResultStruct.barcodeBitmasks
    % hcaClusterConsensusData.clusterResultStruct.barcodeKeys = 
    hcaClusterConsensusData.clusterResultStruct.barcodes = cellfun(@(x) x.stretchedBarcode,barcodeGen(consensusStructs.treeStruct.barMatrix{idx})','UniformOutput',false);
    hcaClusterConsensusData.clusterResultStruct.barcodeBitmasks = cellfun(@(x) x.stretchedrawBitmask,barcodeGen(consensusStructs.treeStruct.barMatrix{idx})','UniformOutput',false);
    hcaClusterConsensusData.clusterResultStruct.flipTFs = logical(consensusStructs.treeStruct.barOrientation{idx}(:,2)-1);
    hcaClusterConsensusData.clusterResultStruct.circShifts = consensusStructs.treeStruct.barOrientation{idx}(:,1)-1;
    hcaClusterConsensusData.clusterResultStruct.barcodeKeys = arrayfun(@(x) num2str(x), find(consensusStructs.treeStruct.barMatrix{idx})', 'UniformOutput', false);
    hcaClusterConsensusData.clusterResultStruct.alignedBarcodes = num2cell(consensusStructs.treeStruct.barcodes{idx},2);
    hcaClusterConsensusData.clusterResultStruct.alignedBarcodeBitmasks = num2cell(~isnan(consensusStructs.treeStruct.barcodes{idx}),2);

    % this is unnecessary..
    hcaClusterConsensusData.details.consensusStruct.clusterResultStructs = {hcaClusterConsensusData.clusterResultStruct};
    hcaClusterConsensusData.details.consensusStruct.formatVersion = '0.1.0';
    hcaClusterConsensusData.details.consensusStruct.timestamp = consensusStruct.time;

    % inputs
    hcaClusterConsensusData.details.consensusStruct.inputs.barcodes =  cellfun(@(x) x.stretchedBarcode,barcodeGen','UniformOutput',false);
    hcaClusterConsensusData.details.consensusStruct.inputs.barcodeBitmasks =   cellfun(@(x) x.stretchedrawBitmask,barcodeGen','UniformOutput',false);
    hcaClusterConsensusData.details.consensusStruct.inputs.barcodeAliases =   cellfun(@(x) x.name,fliplr(kymoStructs)','UniformOutput',false);
    hcaClusterConsensusData.details.consensusStruct.inputs.rawBarcodes =   cellfun(@(x) x.rawBarcode,barcodeGen','UniformOutput',false);
    hcaClusterConsensusData.details.consensusStruct.inputs.rawBgs =   cellfun(@(x) x.rawBg,barcodeGen','UniformOutput',false);
    hcaClusterConsensusData.details.consensusStruct.inputs.clusterThresholdScore = sqrt(length(barcodeGen{1}.stretchedBarcode))*sets.consensus.threshold; %sqrt(N)*clusterthresh

    stretchFactors = cellfun( @(x) x.stretchFactor, barcodeGen);
    for K = 1 : numel(hcaClusterConsensusData.details.consensusStruct.inputs.barcodes)
        hcaClusterConsensusData.details.consensusStruct.inputs.otherBarcodeData{K}.stretchFactor = stretchFactors(K);
        hcaClusterConsensusData.details.consensusStruct.inputs.otherBarcodeData{K}.nmPerPx_original = sets.bitmasking.prestretchPixelWidth_nm;
        hcaClusterConsensusData.details.consensusStruct.inputs.otherBarcodeData{K}.nmPerPx_stretched = sets.bitmasking.prestretchPixelWidth_nm*stretchFactors(K);
        hcaClusterConsensusData.details.consensusStruct.inputs.otherBarcodeData{K}.bpsPerPx_original = NaN;
        hcaClusterConsensusData.details.consensusStruct.inputs.otherBarcodeData{K}.bpsPerPx_stretched = NaN;
    end

    % hcaClusterConsensusData.clusterKey = regexprep(hcaClusterConsensusData.clusterKey,' +',',');

    % key lists. add all barcodes to the beginning of the keylist
    hcaClusterConsensusData.details.consensusStruct.keyList = [arrayfun(@(x)  regexprep(num2str(x),' +',','), find(consensusStructs.treeStruct.barMatrix{idx})', 'UniformOutput', false);...
        cellfun(@(x) regexprep(num2str(cell2mat(x)),' +',',') ,consensusStructs.treeStruct.clusteredBar,'UniformOutput',false)'];
    hcaClusterConsensusData.details.consensusStruct.keyListSimplified = hcaClusterConsensusData.details.consensusStruct.keyList ;
    hcaClusterConsensusData.details.consensusStruct.finalConsensusKey = hcaClusterConsensusData.clusterKey;
    % TODO: add more keys, if there are some unmerged clusters for given selected CC
    % threshold value, otherwise return a single one
    hcaClusterConsensusData.details.consensusStruct.clusterKeys ={hcaClusterConsensusData.clusterKey };

    hcaClusterConsensusData.details.consensusStruct.clusterAssignmentsMatrix = consensusStructs.treeStruct.barMatrix{idx}';
    % next is clusterConsensusData.details.consensusStruct.barcodeStructsMap
    % for each Key, the value is 
    %        maxWeight: 1
    %     indexWeights: [1×212 logical] // how many have been averaged for each
    %     position
    %            alias: 'P18_170607_OD4_100msExp_9-Edited_molecule_1_kymograph.tif Mol #1'
    %          barcode: [1×212 double] // why we need this here too?
    %          parents: {} % this stores the name of previous 2 parents,
    %          orientation and circular shift
    %        bestScore: NaN // not sure what the best score is?
    %      xcorrAtBest: NaN // nan's for first 44 elements since those are just
    %      initial 
    %      
    % mergin tree
    % idx1 idx2 1-xcorr score

    % create container map with values of each cluster key
    ids = [1:length(barcodeGen)];
    names = arrayfun(@(x) num2str(x),ids, 'UniformOutput',false);
    M = containers.Map(names,ids);
    barcodeStructsMap =containers.Map();
    % barcodeStructsMap = containers.Map(names,ids);

    % aadd all barcodes
    for idx=1:length(ids)
        st.indexWeights = barcodeGen{idx}.stretchedrawBitmask;
        st.maxWeight = max(st.indexWeights);
        st.barcode = barcodeGen{idx}.stretchedBarcode;
        st.parents =  {} ;
        st.bestScore = NaN;
        st.xcorrAtBest = NaN;
        st.alias = barcodeGen{idx}.name;
        barcodeStructsMap(num2str(idx)) = st;
    end

    idd = length(barcodeGen); 
    for idx=1:length(consensusStructs.treeStruct.clusteredBar)
        newelt = regexprep(num2str(sort(cell2mat(consensusStructs.treeStruct.clusteredBar{idx}))),' +',',');
        M(newelt)=idd+1;
        idd = idd+1;
        % structure
        st.indexWeights = sum(~isnan(consensusStructs.treeStruct.barcodes{idx}));
        st.maxWeight = max(st.indexWeights);
        st.barcode = nanmean(consensusStructs.treeStruct.barcodes{idx});
        % parents, TF and circshift is not compatible with the old code, since
        % we circularly shift and rotate each molecule, rather than update an
        % average
        st.parents =  {{regexprep(num2str(consensusStructs.treeStruct.clusteredBar{idx}{1}),' +',','),[0],[0]},{regexprep(num2str(consensusStructs.treeStruct.clusteredBar{idx}{2}),' +',','),[0],[0]} } ;
        st.bestScore = NaN;
        st.alias ='';

        st.xcorrAtBest = NaN;

        barcodeStructsMap(newelt) = st;
    end
    hcaClusterConsensusData.details.consensusStruct.barcodeStructsMap = barcodeStructsMap;

    % merging tree
    w1= cellfun(@(x) M(regexprep(num2str(sort(x{1})),' +',',')),consensusStructs.treeStruct.clusteredBar)';
    w2 = cellfun(@(x) M(regexprep(num2str(sort(x{2})),' +',',')),consensusStructs.treeStruct.clusteredBar)';
    hcaClusterConsensusData.details.consensusStruct.consensusMergingTree = [w1 w2  1-consensusStructs.treeStruct.maxCorCoef'];


    hcaClusterConsensusData.details.consensusStruct.timestamp = consensusStruct.time;
    hcaClusterConsensusData.details.consensusStruct.formatVersion = '0.1.0'; % is repeating

    % cache saves all the xcorrs.
    hcaClusterConsensusData.details.consensusStruct.cache = containers.Map();
    clusterConsensusData = hcaClusterConsensusData;
    save(fullfile(sets.output.matDirpath,'cluster.mat'),'-v7.3','clusterConsensusData')


end

