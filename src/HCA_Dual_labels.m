function [hcaStruct] = HCA_Dual_labels(sets, hcaStruct)
%% HCA_Dual_labels
% Used for comparing fagments of human chromosome to chromosomes using
% CB (competitive binding) and flourescent markers (dots).
%
%     Args:
%         sets (struct): Input settings to the method
%         hcaStruct (struct): Input structure, if non-empty, load
%         result structure instead of computing everything
%
%     Returns:
%         hcaStruct: Return structure
%
%     Example:
%         This is an example: run [hcaStruct] = HCA_Dual_labels(sets)
%

if nargin < 1 % if settings were not provided
  %% import settings
  import CBT.Hca.Import.import_hca_settings;
  sets = import_hca_settings('hca_settings_dual_label.txt');
  assignin('base','sets', sets)
end

if nargin < 2
  %% Ask for barcodes and import them
  import DL.Hca.import_single_timeframe_barcodes
  barcodeGenDense = import_single_timeframe_barcodes(sets);
  if size(barcodeGenDense, 1) > 1
    barcodeGen = barcodeGenDense;
  else
    import DL.Hca.import_sparse_barcodes
    barcodeGenSparse = import_sparse_barcodes(barcodeGenDense, sets);
    import DL.Hca.merge_barcodegens
    barcodeGen = merge_barcodegens(barcodeGenDense, barcodeGenSparse);
  end
  assignin('base','barcodeGen', barcodeGen)
  
  
  %% now user theories. They could already be in txt files (if generated
  % (with HCA 4.0.0), but we should keep support for older theory files too
  sets.theoryFile=[];
  sets.theoryFileFold = [];
  
  % get user theory
  import DL.Hca.get_user_theory_dl;
  theoryStruct = get_user_theory_dl(sets);
  assignin('base','theoryStruct', theoryStruct)
  
  %   import CBT.Hca.Settings.get_user_theory;
  %   theoryStruct = get_user_theory(sets);
  %   tmptheoryStruct = get_user_theory(sets);
  %   theoryStruct.filename2 = tmptheoryStruct.filename;
  %   assignin('base','theoryStruct', theoryStruct)
  
  [molTableFileName, dirpath] = uigetfile('*.*', 'Select file containing molecule data-table');
  externalAlignmentStruct = struct();
  externalAlignmentStruct.idx = 1;
  if not(isequal(dirpath, 0))
    molDataTable = readtable(fullfile(dirpath, molTableFileName));
    externalAlignmentStruct.molIds = cellfun(@(x) str2double(x.name(1:end-4)), barcodeGen(1,:));
    [~, tableIds] = ismember(externalAlignmentStruct.molIds, molDataTable.MoleculeID);
    externalAlignmentStruct.Z = -norminv(10.^-molDataTable.Confidence(tableIds));
    refStartPosBp = min(molDataTable.RefStartPos(tableIds), molDataTable.RefEndPos(tableIds));
    qryFirstDotPos = min(molDataTable.QryStartPos(tableIds), molDataTable.QryLen(tableIds) - molDataTable.QryStartPos(tableIds));
    externalAlignmentStruct.firstDotPos = qryFirstDotPos./molDataTable.QryLen(tableIds).*molDataTable.Length(tableIds)/375;
    externalAlignmentStruct.pos = max(1, ceil(refStartPosBp/sets.pvalue.pixelWidth_nm*sets.pvalue.nmbp-externalAlignmentStruct.firstDotPos));
    externalAlignmentStruct.or = strcmp(molDataTable.Orientation(tableIds), '-')+1;
    if length(theoryStruct) > 1
      list = cellfun(@(x) x.name, theoryStruct, 'un', 0);
      externalAlignmentStruct.idx = listdlg('PromptString', ...
        {'Select which theory the molecule' 'positions are referring to:'}, ...
        'ListString', list, 'SelectionMode', 'Single');
    end
  end
  assignin('base','externalAlignmentStruct', externalAlignmentStruct)
  
  %% compare theory to experiment
  rezMax = cell(1,length(theoryStruct));
  import DL.Hca.on_compare_dual_label;
  for barNr = 4%1:length(theoryStruct)
    disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));
    if not(isempty(fieldnames(externalAlignmentStruct))) && externalAlignmentStruct.idx == barNr
      rezMax{barNr} = on_compare_dual_label( ...
        barcodeGen, ...
        theoryStruct{barNr}, ...
        20, ...
        sets, ...
        externalAlignmentStruct);
    else
      rezMax{barNr} = on_compare_dual_label( ...
        barcodeGen, ...
        theoryStruct{barNr}, ...
        20, ...
        sets);
    end
  end
  assignin('base','rezMax', rezMax)
  
  import DL.Hca.combine_theory_results_dl;
  comparisonStruct = combine_theory_results_dl( ...
    theoryStruct, ...
    rezMax);
  assignin('base','comparisonStruct', comparisonStruct)
  
  import DL.Hca.get_display_results_dl;
  get_display_results_dl( ...
    barcodeGen, ...
    comparisonStruct, ...
    theoryStruct, ...
    sets, ...
    externalAlignmentStruct);
  
  
else
  % run load results function or something similar..
end

end
