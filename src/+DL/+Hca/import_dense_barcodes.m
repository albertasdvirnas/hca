function denseStruct = import_dense_barcodes(sets)
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
  importStruct{i}.name = itemsToImportFilenames{i};
  importStruct{i}.alignedKymo = im2double(imread(fullfile(itemsToImportFolders{i}, ...
    itemsToImportFilenames{i})));
  importStruct{i}.leftEdgeIdxs = 1;
  importStruct{i}.rightEdgeIdxs = length(importStruct{i}.alignedKymo);
end

% generate barcodes
import CBT.Hca.Core.gen_barcodes;
denseStruct = CBT.Hca.Core.gen_barcodes(importStruct, sets);

import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode
for i=1:length(denseStruct)
  denseStruct{i}.rawBarcode = filter_barcode( ...
    denseStruct{i}.rawBarcode, ...
    sets.filterSettings);
end
end

