function sparseStruct = import_sparse_barcodes(barcodeGen, sets)
% loads figure window
import Fancy.UI.Templates.create_figure_window;
import DL.Hca.create_import_tab;

answer = questdlg('Is the second type of barcode an intensity profile or a list of dot positions?', ...
  'Choose label type', ...
  'Intensity profile', ...
  'Dots', ...
  'Dots');
if strcmp(answer, 'Intensity profile')
  import DL.Hca.import_dense_barcodes
  sparseStruct = import_dense_barcodes(sets);
  return
end

cache = containers.Map();
while true
  [hMenuParent, tsHCA] = create_figure_window( ...
    "Import dots", ...
    'Dual-label HCA');
  cache = create_import_tab(hMenuParent,tsHCA,'dots');
  uiwait(gcf);
  delete(hMenuParent);
  if not(size(cache('selectedItems'),1) == length(barcodeGen))
    f = msgbox('The number of dot-maps did not match the number of barcodes');
    uiwait(f);
  else
    break
  end
end
itemsToImport = cache('selectedItems');
namesToSort = itemsToImport(:,1);
[~, sid] = sort(namesToSort);
itemsToImportFilenames = itemsToImport(sid,1);
itemsToImportFolders = itemsToImport(sid,2);

import Microscopy.Simulate.Core.apply_point_spread_function
% Temp struct for gen_barcodes
importStruct = cell(1,length(itemsToImportFilenames));
for i=1:length(itemsToImportFilenames)
  thisDots = importdata(fullfile(itemsToImportFolders{i}, ...
    itemsToImportFilenames{i}));
  dotbar = zeros(1, length(barcodeGen{i}.rawBarcode));
  dotbar(ceil(thisDots)) = sets.bitmasking.prestretchPixelWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
  importStruct{i}.name = itemsToImportFilenames{i};
  importStruct{i}.alignedKymo = apply_point_spread_function( ...
    dotbar, ...
    sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm, ...
    1);
  importStruct{i}.leftEdgeIdxs = 1;
  importStruct{i}.rightEdgeIdxs = length(importStruct{i}.alignedKymo);
end

% generate barcodes
import CBT.Hca.Core.gen_barcodes;
sparseStruct = CBT.Hca.Core.gen_barcodes(importStruct, sets);
end

