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
  sets = import_hca_settings('hca_settings.txt');
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
  
  %% compare theory to experiment
  import CBT.Hca.Core.Comparison.compare_distance;
  rezMax = compare_distance( ...
    barcodeGen, ...
    theoryStruct, ...
    sets);
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
    sets);
  
  
else
  % run load results function or something similar..
end

end
