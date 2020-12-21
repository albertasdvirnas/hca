function sparseStruct = import_sparse_barcodes(barcodeGen, sets)
% import sparse barcodes

% If
if ~sets.kymosets.askforkymos
  try
    sets.kymosets.kymoFile = sets.kymosets.sparseFile;
    if isequal(sets.kymosets.sparseMap, 1)
      import DL.Hca.import_single_timeframe_barcodes
      sparseStruct = import_single_timeframe_barcodes(sets);
      return
    end
  catch
    sets.kymosets.askforkymos = 1;
  end
end

if sets.kymosets.askforkymos
  % loads figure window
  import Fancy.UI.Templates.create_figure_window;
  import DL.Hca.create_import_tab;
  answer = questdlg('What is the second type of barcode?', ...
    'Choose label type', ...
    'Intensity profile', ...
    'Dots', ...
    'Dots');
  if strcmp(answer, 'Intensity profile')
    import DL.Hca.import_single_timeframe_barcodes
    sparseStruct = import_single_timeframe_barcodes(sets);
    return
  end
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
  dotbar(min(ceil(thisDots), length(dotbar))) = sets.bitmasking.prestretchPixelWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
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

