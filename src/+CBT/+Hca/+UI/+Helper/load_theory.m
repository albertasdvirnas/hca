function [ theoryStruct ] = load_theory( sets )
    % load_theory
    % Loads theory names from path
    %     Args:
    %         barcodeGen: barcode structure
    %         sets: Input settings to the method
    % 
    %     Returns:
    %         consensusStructs: Consensus tree structure
    % 
    %     Example:
    %          
    
    % todo: make it recognize that the file list is txt files,
    % and then just save these names to txt struct
    
    % In case theory is in a single file, read of it's contents into a
    % folder and then save the names of the theories
    
   % theoryStruct = cell(1,length(barcodeData.hcaSessionStruct.theoryBarcodes));
    theoryStruct = {};

    for idx=1:length(sets.theoryFileFold)
        [~,part,st] = fileparts(sets.theoryFile{idx});
        if isequal(st,'.txt')
            % file already exists, just load the info from header
            names = fullfile(sets.theoryFileFold{idx},sets.theoryFile{idx});
            parts = strsplit(names,'_');
            theoryStructN.filename = names;
            theoryStructN.name = part;
            theoryStructN.length = round(str2num(parts{end-5}));
            theoryStructN.meanBpExt_nm  = str2num(parts{end-4});
            theoryStructN.pixelWidth_nm =  str2num(parts{end-3});
            theoryStructN.psfSigmaWidth_nm = str2num(parts{end-2});
            theoryStructN.isLinearTF = str2num(parts{end-1});
            theoryStruct = [theoryStruct theoryStructN];
        else
            %         addpath(genpath(sets.theoryFileFold{idx}));
            barcodeData = load(fullfile(sets.theoryFileFold{idx},sets.theoryFile{idx}));

            try % there can be different kind of input data, so use the appropriate one
                bars = barcodeData.hcaSessionStruct.theoryGen.theoryBarcodes;
                names = barcodeData.hcaSessionStruct.theoryGen.theoryNames;
                setsB = barcodeData.hcaSessionStruct.theoryGen.sets;
            catch
                try
                    bars = barcodeData.hcaSessionStruct.theoryBarcodes;
                    names = barcodeData.hcaSessionStruct.theoryNames;
                    setsB = barcodeData.hcaSessionStruct.sets;
                catch
                    try
                        bars = barcodeData.theoryGen.theoryBarcodes;
                        names = barcodeData.theoryGen.theoryNames;
                        setsB = barcodeData.theoryGen.sets;        
                    catch
                        bars = barcodeData.theoryBarcodes;
                        names = barcodeData.theoryNames;
                        setsB = barcodeData.sets; 
                    end
                end 
            end
        
            theoryStructN = cell(1,length(bars));
            theoryFold = sets.theoryFileFold{idx};
            precision = sets.theory.precision;
            % these settings from old file // if all files are txts, read this
            % from filename
            meanbpnm=setsB.meanBpExt_nm;
            pixelWidth_nm = setsB.pixelWidth_nm;
            psfSigmaWidth_nm = setsB.psfSigmaWidth_nm;
            try  % try to see if setting of linear theory function is included
               isLinearTF  =  setsB.psfSigmaWidth_nm;
            catch
               isLinearTF  =  0;
            end

            parfor i=1:length(bars)
                fname = fullfile(theoryFold, strcat(num2str(i),'.txt'));
                fileID = fopen(fname,'w');
                % choose the precision from settings
                fprintf(fileID,strcat(['%2.' num2str(precision) 'f ']),bars{i});
                fclose(fileID);
                theoryStructN{i}.filename = fname;
                theoryStructN{i}.name = names{i};
                theoryStructN{i}.length =  length(bars{i});
                theoryStructN{i}.meanBpExt_nm = meanbpnm;
                theoryStructN{i}.pixelWidth_nm = pixelWidth_nm;
                theoryStructN{i}.psfSigmaWidth_nm = psfSigmaWidth_nm;
                theoryStructN{i}.isLinearTF = isLinearTF;
            end
            theoryStruct = [theoryStruct theoryStructN];
        end
    %
    end
end

