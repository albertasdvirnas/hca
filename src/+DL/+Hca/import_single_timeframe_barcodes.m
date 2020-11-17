function barcodeGen = import_single_timeframe_barcodes(sets)
    % Args: sets: settings file
    % Returns: barcodeGen - barcode structure file


    % allow to pass a file with single timeframe barcodes instead of 
    % loading a GUI
    if ~sets.kymosets.askforkymos 
        try 
            fid = fopen(sets.kymosets.kymoFile); 
            fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
            for i=1:length(fastaNames{1})
                [FILEPATH,NAME,EXT] = fileparts(fastaNames{1}{i});
                itemsToImportFilenames{i} = strcat(NAME,EXT);
                itemsToImportFolders{i} = FILEPATH;
            end
        catch
            sets.kymosets.askforkymos   = 1;
        end
    end

    
    if sets.kymosets.askforkymos 
        % loads figure window
        import Fancy.UI.Templates.create_figure_window;
        import DL.Hca.create_import_tab;

        cache = containers.Map();
        while true
          [hMenuParent, ...
            tsHCA] = create_figure_window( ...
            "Import barcodes", ...
            'Dual-label HCA');
          cache = create_import_tab(hMenuParent,tsHCA,'barcodes',cache);
          uiwait(gcf);
          delete(hMenuParent);
          break
        end
        itemsToImport = cache('selectedItems');
        namesToSort = itemsToImport(:,1);
        [~, sid] = sort(namesToSort);
        itemsToImportFilenames = itemsToImport(sid,1);
        itemsToImportFolders = itemsToImport(sid,2);
    end

    % Temp struct for gen_barcodes
    importStruct = cell(1,length(itemsToImportFilenames));
    for i=1:length(itemsToImportFilenames)
        filePath = fullfile(itemsToImportFolders{i}, itemsToImportFilenames{i});
        [~, ~, fileExt] = fileparts(filePath);
        switch fileExt
            case '.mat'
                tmpMatStruct = load(filePath);
                tmpMatFields = fields(tmpMatStruct);
                if length(itemsToImportFilenames) > 1
                    importStruct{i} = subsref(tmpMatStruct, substruct('.', tmpMatFields{1}));
                else
                    barcodeGen = subsref(tmpMatStruct, substruct('.', tmpMatFields{1}));
                    return;
                end
            case '.txt'
                importStruct{i}.alignedKymo = importdata(filePath);
            otherwise
                importStruct{i}.alignedKymo = im2double(imread(filePath));
        end
        importStruct{i}.name = itemsToImportFilenames{i};
        importStruct{i}.leftEdgeIdxs = 1;
        importStruct{i}.rightEdgeIdxs = length(importStruct{i}.alignedKymo);

    end

    % generate barcodes
    import CBT.Hca.Core.gen_barcodes;
    barcodeGen = CBT.Hca.Core.gen_barcodes(importStruct, sets);

    import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode
    for i=1:length(barcodeGen)
      barcodeGen{i}.rawBarcode = filter_barcode( ...
        barcodeGen{i}.rawBarcode, ...
        sets.filterSettings);
    end
    
end

