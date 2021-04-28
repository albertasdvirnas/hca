function barcodeGen = merge_barcodegens(barcodeGenDense, barcodeGenSparse)
for i=1:length(barcodeGenDense)
  if isempty(barcodeGenDense{i}) || isempty(barcodeGenSparse{i})
    throw(MException('barcodeGenMerge:empty', ...
      'Barcode pair member is empty, aborting barcode import.'));
  end
  if not(strcmp(barcodeGenDense{i}.name(1:end-4), ...
      barcodeGenSparse{i}.name(1:end-4)))
    warning('Difference found in barcode pair names.')
  end
%   if not(length(barcodeGenSparse{i}.rawBarcode) ...
%       == length(barcodeGenDense{i}.rawBarcode))
%     throw(MException('barcodeGenMerge:lengthDiscrepancy', ...
%       'Barcode length discrepancy detected, aborting barcode import.'));
%   end
barcodeGen = [barcodeGenDense; barcodeGenSparse];
end

