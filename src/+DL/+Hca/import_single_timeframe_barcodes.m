function barcodeGen = import_single_timeframe_barcodes(sets)
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
        return
      end
    otherwise
      importStruct{i}.name = itemsToImportFilenames{i};
      importStruct{i}.alignedKymo = im2double(imread(filePath));
      importStruct{i}.leftEdgeIdxs = 1;
      importStruct{i}.rightEdgeIdxs = length(importStruct{i}.alignedKymo);
  end
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

