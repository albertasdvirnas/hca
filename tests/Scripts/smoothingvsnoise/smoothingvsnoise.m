%
k=300;
numBarcodes = 10;



% generate a hundred series
rV =  normrnd(0,1,numBarcodes,k);

par = 20;
% what happens if they A and B have different filtering?
for ix=1:size(rV,1)
    rV(ix,:) = imgaussfilt(rV(ix,:),par);
end

% now compare all of the to each other to find the best CC

% dist = MASS_PCC([rV(2,:) rV(2,1:end-1)],rV(1,:),k);

maxcoef = zeros(numBarcodes-1,numBarcodes-1);
or = zeros(numBarcodes-1,numBarcodes-1);
pos = zeros(numBarcodes-1,numBarcodes-1);



% run from the first barcode to the almost last one
for barcodeIdxA = 1:numBarcodes-1
    barcodeA =  rV(barcodeIdxA,:);
    for barcodeIdxB = barcodeIdxA+1:numBarcodes
        barcodeB = [rV(barcodeIdxB,:) rV(barcodeIdxB,1:end-1)];
        % compute correlatiton coefficients
        xcorrs = MASS_PCC(barcodeB,barcodeA,k);

        % compute maximum, position and orientation for the best
        % position
        [f,s] = max(xcorrs);          
        [ b, ix ] = sort( f(:), 'descend' );
        maxcoef(barcodeIdxA,barcodeIdxB-1) = b(1);
        or(barcodeIdxA,barcodeIdxB-1) = s(ix(1));
        pos(barcodeIdxA,barcodeIdxB-1) = ix(1);
    end     
end

colvect= maxcoef(find(~tril(ones(size(maxcoef)))));

    max(colvect)
    mean(colvect)