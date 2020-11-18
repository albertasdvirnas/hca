function [t,matFilepathShort,theoryStruct, sets,theoryGen] = HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

% This function computes theory barcode for a sequence
% load default settings

% timestamp for the results
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

if nargin < 3 
    setsFile = '';
end

% import settings
import CBT.Hca.Import.import_settings;
[sets] = import_settings(setsFile);

if nargin < 1
    sets.fastas = 'theories_parallel.txt';
else
    if nargin >= 5
        data = cellfun(@(x,y) fullfile(x,y),theoryfold,theoryfile,'UniformOutput',false);
        sets.theoryNames = theoryfile;
    else    
        listing = dir(theory_names);
        data = arrayfun(@(x) fullfile(x.folder,x.name),listing,'UniformOutput',false);
        
    end
    fd = fopen('theories_parallel.txt','w');
    for i=1:length(data)
        fprintf(fd,'%s\n',data{i});
    end
    fclose(fd);
    sets.fastas = 'theories_parallel.txt';
end

% load default settings
import CBT.Hca.Settings.get_user_theory_sets;
sets = get_user_theory_sets(sets);


if ~sets.skipBarcodeGenSettings
    import CBT.Hca.Settings.get_theory_sets;
    sets.theoryGen = get_theory_sets(sets.theoryGen); %
end

tic;
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
theoryBarcodes = cell(1,length(sets.theoryNames));
theoryNames = cell(1,length(sets.theoryNames));
theoryIdx= cell(1,length(sets.theoryNames));
bpNm= cell(1,length(sets.theoryNames));
tempNames = cell(1,length(sets.theoryNames));

theories = sets.theoryFold;
theorynames = sets.theoryNames;

if nargin < 2
    meanBpExt_nm = sets.theoryGen.meanBpExt_nm;
end
sets.theoryGen.meanBpExt_nm = meanBpExt_nm;
bpNmV = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;

meanBpExt_nm = sets.theoryGen.meanBpExt_nm;
pixelWidth_nm = sets.theoryGen.pixelWidth_nm;
psfSigmaWidth_nm = sets.theoryGen.psfSigmaWidth_nm;
linear = sets.theoryGen.isLinearTF;
resultsDir = sets.resultsDir;
% loop over theory file folder
parfor idx = 1:length(sets.theoryNames)

%     addpath(genpath(theories{idx}))
    disp(strcat(['loaded theory sequence ' theorynames{idx}] ));

    % new way to generate theory, check theory_test.m to check how it works
    import CBT.Hca.Core.Theory.compute_theory_barcode;
    [theorySeq, header,bitmask] = compute_theory_barcode(fullfile(theories{idx},theorynames{idx}),sets);

	theoryBarcodes{idx} = theorySeq;
    theoryBitmasks{idx} = bitmask;

    theoryNames{idx} = header;
    theoryIdx{idx} = idx;
    bpNm{idx} = bpNmV;
    
    
        
    if sets.savetxts && ~isempty(bitmask)
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        tempNames{idx} = strcat(['theory_' C{1} '_' num2str(length(bitmask)) '_' num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_bitmask.txt']);
%         matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile(resultsDir, timestamp, tempNames{idx});
        fd = fopen(matFilepath,'w');
        fprintf(fd, strcat([' %5.5f']), bitmask);
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
theoryGen.theoryBarcodes = theoryBarcodes;
theoryGen.theoryBitmasks = theoryBitmasks;

theoryGen.theoryNames = theoryNames;
theoryGen.theoryIdx = theoryIdx;
theoryGen.bpNm = bpNm;

theoryGen.sets = sets.theoryGen;

matFilename = strcat(['theoryStruct_' num2str(meanBpExt_nm) '_' sprintf('%s_%s', timestamp) 'session.mat']);
matFilepath = fullfile(sets.resultsDir, matFilename);

save(matFilepath, 'theoryGen');


matTpathShort = fullfile(sets.resultsDir, strcat(['theories_mat_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matTpathShort,'w');    
fprintf(fd, '%s \n',matFilepath);
fclose(fd);

% print out the results
fprintf('Saved theory fasta names ''%s'' to ''%s''\n', matFastapathShort, matFilepath);
% fprintf('Saved theory txts names ''%s'' to ''%s''\n', matFilepathShort, matFilepath);
% fprintf('Saved theory struct data ''%s'' to ''%s''\n', matFilename, matFilepath);
fprintf('Saved theory mat filename ''%s'' to ''%s''\n', matFilename, matTpathShort);
t=toc;

addpath(genpath(sets.resultsDir));

sets.theories = matFilepathShort;
sets.theory.askfortheory = 0;
sets.theory.askfornmbp = 0;
sets.theory.nmbp = sets.theoryGen.meanBpExt_nm;
sets.theory.precision = 5;
sets.theory.askfortheory = 0;
% Now compare vs theory
% get user theory
import CBT.Hca.Settings.get_user_theory;
[theoryStruct, sets] = get_user_theory(sets);
    
end

