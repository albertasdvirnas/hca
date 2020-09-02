function [barcodeGenC] = gen_subfragments(barcodeGen,numberFragments)

    jj = 1;
    barcodeGenC = cell(1,numberFragments*length(barcodeGen));
    for i=1:length(barcodeGen)
        lengthSubfragment = round(sum(barcodeGen{i}.rawBitmask)/numberFragments);

        st = find(barcodeGen{i}.rawBitmask,1,'first');

        for j=1:numberFragments
            barcodeGenC{jj} = barcodeGen{i};
            barcodeGenC{jj}.rawBarcode = barcodeGen{i}.rawBarcode(st+lengthSubfragment*(j-1):st+lengthSubfragment+lengthSubfragment*(j-1)-1);
            barcodeGenC{jj}.rawBitmask = ones(1,length(barcodeGenC{jj}.rawBarcode ));
            parts = strsplit(barcodeGen{i}.name,'_');
            parts{2} = num2str(str2num(parts{2})+st+lengthSubfragment*(j-1)-1);
            parts{3} = num2str(lengthSubfragment);
            barcodeGenC{jj}.name = strjoin([parts,num2str(jj)],'_');
%             barcodeGenC{jj}.name = strcat([num2str(j) '_' barcodeGenC{i}.name]);
            jj = jj+1;
        end

    end
            
end

