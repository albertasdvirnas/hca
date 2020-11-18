% This function computes theory barcode for a sequence
% load default settings

% timestamp for the results
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% import settings
import CBT.Hca.Import.import_settings;
[sets] = import_settings('theory_settings.txt');

% load default settings
import CBT.Hca.Settings.get_user_theory_sets;
sets = get_user_theory_sets(sets);


if ~sets.skipBarcodeGenSettings
    import CBT.Hca.Settings.get_theory_sets;
    sets.theoryGen = get_theory_sets(sets.theoryGen); %
end

theories = sets.theoryFold;
theorynames = sets.theoryNames;


% make theoryData folder
mkdir(sets.resultsDir);
mkdir(sets.resultsDir,timestamp);

% file where the names of theories generated in this session will be saved
matFilepathShort = fullfile(sets.resultsDir, strcat(['theories_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matFilepathShort,'w');    
fclose(fd);

matFastapathShort = fullfile(sets.resultsDir, strcat(['fastas_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matFastapathShort,'w');    
fclose(fd);

% compute free concentrations
import CBT.Hca.Core.Theory.compute_free_conc;
sets = compute_free_conc(sets);

theoryGen = struct();
tempNames = cell(1,length(sets.theoryNames));

% sets.theoryGen.meanBpExt_nm = meanBpExt_nm;
bpNmV = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;

meanBpExt_nm = sets.theoryGen.meanBpExt_nm;
pixelWidth_nm = sets.theoryGen.pixelWidth_nm;
psfSigmaWidth_nm = sets.theoryGen.psfSigmaWidth_nm;
linear = sets.theoryGen.isLinearTF;
resultsDir = sets.resultsDir;

% loop over theory file folder
for idx = 1:length(sets.theoryNames)

    addpath(genpath(sets.theoryFold{idx}))
    disp(strcat(['loaded theory sequence ' sets.theoryNames{idx}] ));

    % new way to generate theory, check theory_test.m to check how it works
    import CBT.Hca.Core.Theory.compute_theory_barcode;
    [theorySeq, header,bitmask] = compute_theory_barcode(sets.theoryNames{idx},sets);

	theoryGen.theoryBarcodes{idx} = theorySeq;
    theoryGen.theoryBitmasks{idx} = bitmask;

    theoryGen.theoryNames{idx} = header;
    theoryGen.theoryIdx{idx} = idx;
    theoryGen.bpNm{idx} = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;
    
     
        
    if sets.savetxts && ~isempty(bitmask)
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        tempNames{idx} = strcat(['theory_' C{1} '_' num2str(length(bitmask)) '_' num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_bitmask.txt']);
%         matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile(resultsDir, timestamp, tempNames{idx});
        fd = fopen(matFilepath,'w');
        fprintf(fd,' %5d', bitmask);
        fclose(fd);
    end
    
    if sets.savetxts && ~isempty(theorySeq)
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        tempNames{idx} = strcat(['theory_' C{1} '_' num2str(length(theorySeq)) '_' num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_barcode.txt']);
%         matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile(resultsDir, timestamp, tempNames{idx});
        fd = fopen(matFilepath,'w');
        fprintf(fd, strcat([' %5.' num2str(sets.theoryGen.precision) 'f ']), theorySeq);
        fclose(fd);
    end
end


for idx = 1:length(sets.theoryNames)
    if sets.savetxts
        fd = fopen(matFilepathShort,'a'); fprintf(fd, '%s \n',fullfile(resultsDir,timestamp, tempNames{idx})); fclose(fd);
        fd = fopen(matFastapathShort,'a'); fprintf(fd, '%s \n',fullfile(theories{idx},theorynames{idx})); fclose(fd);
    end
end

% save sets
theoryGen.sets = sets.theoryGen;

matFilename = strcat(['theoryStruct_' sprintf('%s_%s', timestamp) 'session.mat']);
matFilepath = fullfile(sets.resultsDir, matFilename);

save(matFilepath, 'theoryGen');


matTpathShort = fullfile(sets.resultsDir, strcat(['theories_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matTpathShort,'w');    
fprintf(fd, '%s \n',matFilepath)
fclose(fd);

% print out the results
fprintf('Saved theory fasta names ''%s'' to ''%s''\n', matFastapathShort, matFilepath);
fprintf('Saved theory txts names ''%s'' to ''%s''\n', matFilepathShort, matFilepath);
fprintf('Saved theory struct data ''%s'' to ''%s''\n', matFilename, matFilepath);
fprintf('Saved theory mat filename ''%s'' to ''%s''\n', matFilename, matTpathShort);

    
