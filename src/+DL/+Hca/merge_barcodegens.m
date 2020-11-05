function barcodeGen = merge_barcodegens(barcodeGenDense, barcodeGenSparse)
for i=1:length(barcodeGenDense)
  if not(strcmp(barcodeGenDense{i}.name, barcodeGenSparse{i}.name))
    warning('Difference found in barcode pair names.')
  end
  if not(length(barcodeGenSparse{i}.rawBarcode) ...
      == length(barcodeGenDense{i}.rawBarcode))
    error('Barcode length anomaly detected, aborting barcode import.');
  end
  barcodeGenDense{i}.rawDotBar = barcodeGenSparse{i}.rawBarcode;
end
barcodeGen = barcodeGenDense;

