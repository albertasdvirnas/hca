function [hcaStruct] = HCA_Dual_labels(sets, hcaStruct)
% HCA_Dual_labels
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
  % import settings
  import CBT.Hca.Import.import_hca_settings;
  sets = import_hca_settings('hca_settings.txt');
  assignin('base','sets', sets)
end

if nargin < 2
  % Ask for barcodes and import them  
  import DL.Hca.import_single_timeframe_barcodes
  barcodeGenDense = import_single_timeframe_barcodes(sets);
  if size(barcodeGenDense, 1) > 1
    barcodeGen = barcodeGenDense;
  else
    sets.kymosets.kymoFile = sets.kymosets.sparseFile;
    import DL.Hca.import_sparse_barcodes
    barcodeGenSparse = import_sparse_barcodes(barcodeGenDense, sets);
    import DL.Hca.merge_barcodegens
    barcodeGen = merge_barcodegens(barcodeGenDense, barcodeGenSparse);
  end
  
  
  %% now user theories. They could already be in txt files (if generated
  % (with HCA 4.0.0), but we should keep support for older theory files too
  sets.theoryFile=[];
  sets.theoryFileFold = [];
  
  % get user theory
  import CBT.Hca.Settings.get_user_theory;
  [theoryStruct, ~] = get_user_theory(sets);
  [tmpTheoryStruct, ~] = get_user_theory(sets);
  for i=1:length(theoryStruct)
    theoryStruct{i}.filename2 = tmpTheoryStruct{i}.filename;
  end
  assignin('base','theoryStruct', theoryStruct)
  
  % compare theory to experiment
  import CBT.Hca.Core.Comparison.compare_distance;
  [rezMax, bestBarStretch, bestLength] = compare_distance( ...
    barcodeGen, ...
    theoryStruct, ...
    sets);
  
  import CBT.Hca.Core.Comparison.combine_theory_results;
  comparisonStruct = combine_theory_results( ...
    theoryStruct, ...
    rezMax, ...
    bestBarStretch, ...
    bestLength);
  assignin('base','comparisonStruct', comparisonStruct)
  
  import DL.Hca.get_display_results_dl;
  get_display_results_dl( ...
    barcodeGen, ...
    struct(), ...
    comparisonStruct, ...
    theoryStruct, ...
    sets);
  
%   import CBT.Hca.Core.additional_computations
%   additional_computations( ...
%     barcodeGen, ...
%     struct(), ...
%     comparisonStruct, ...
%     theoryStruct, ...
%     comparisonStructAll, ...
%     sets);
  
  
  
else
  % run load results function or something similar..
end

end
