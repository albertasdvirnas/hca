function [ rezMaxM,bestBarStretch,bestLength ] = on_compare(barcodeGen,theoryStruct,comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask)
    % on_compare_theory_to_exp
    % Compares experiments to single theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theoryStruct: theory structure
    % 
    %     Returns:
    %         comparisonStructure: comparison structure

    % create the rezult structure
%     comparisonStructure = cell(1, length(barcodeGen));
    %ccM = cell(length(rawBarcodes),1);

    switch comparisonMethod
        case 'spearman'
            import SignalRegistration.unmasked_spearman_corr;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr(x',y',w);
        case 'spearman2'
            import SignalRegistration.unmasked_spearman_corr2;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr2(x',y',w);
        case 'spearman3'
            import SignalRegistration.unmasked_spearman_corr3;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr3(x',y',w);

        case 'unmasked_pcc_corr'
            import SignalRegistration.unmasked_pcc_corr;
            comparisonFun = @(x,y,z,w,u) unmasked_pcc_corr(x,y,z);
%             import CBT.Hca.UI.Helper.get_best_parameters;

%             parameterfun = @(x)
        case 'mass_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);
        case 'dtw'
            
            import SignalRegistration.ucr_dtw_score;
            comparisonFun = @(X1,X2, bitX1, bitX2, w) ucr_dtw_score(X1,X2, bitX1);

%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets);
        case 'mp'
            import mp.mp_dist_stomp_with_masks;
            comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),[],theoryStruct.isLinearTF);

        case 'mpAll' % mp for all subfragments separately
            import mp.mp_dist_stomp_with_masks_all;
             comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks_all(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))));

%             comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks_all(X1',[X2 X2(1:w-1)]', bitX1', bitX2', w,2^(4+nextpow2(length(X1))));

%             [scoreMatPCC,scorePos,orientation,secondPos, abMatchTemp] = mp_dist_stomp_with_masks(X1,X2, bitX1, bitX2, w, kk,par,islinear)

%             comparisonFun = @(x,y,z,w) unmasked_MP(y,x,z,2^(4+nextpow2(length(x))));
% 
%             error('not yet implemented');
            
        case 'hmm'
            % here we stor the hmm method used for structural variation
            % detection. However, the outputs need to be aligned...
            error('not yet implemented..');

        case 'mpdist'
            error('not yet implemented');
        otherwise
            error('undefined method');
    end
            

    % stretch factors
%     stretchFactors = sets.theory.stretchFactors;
%     try % some versions of sets might not have this number, therefore we add a test here
%         minLength = sets.comparison.minLength;
%     catch
%         minLength = 20;
%     end
    
    % number of pixels around best theory match coefficient
%     try
%         numPixelsAroundBestTheoryMask = sets.bitmasking.numPixels;
%     catch
%         numPixelsAroundBestTheoryMask = 0;
%     end

    
    % load theory barcode txt file. For UCR DTW (c++ code), we only need the name of the
    % file so this can be skipped.
    fileID = fopen(theoryStruct.filename,'r');
    formatSpec = '%f';
    theorBar = transpose(fscanf(fileID,formatSpec));
    fclose(fileID);
    theorBit = ones(1,length(theorBar));
    
    % if we want match to be circular

   	import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode

    rezMaxM = cell(1,length(barcodeGen));
    bestBarStretch = zeros(1,length(barcodeGen));
    bestLength = zeros(1,length(barcodeGen));
    % for all the barcodes run
    for idx=1:length(barcodeGen)
%         idx
        % xcorrMax stores the  maximum coefficients
        xcorrMax = zeros(1,length(stretchFactors));
        
        % rezMaz stores the results for one barcode
        rezMax = cell(1,length(stretchFactors));
       
        try
        % barTested barcode to be tested
        barTested = barcodeGen{idx}.rawBarcode;
        catch
        barTested = barcodeGen{idx}.barcode;
        end
        
        % in case barcode should be filtered
%         barTested = filter_barcode(barTested, sets.filterSettings);

        % length of this barcode
        lenBarTested = length(barTested);
        
        % barBitmask - bitmask of this barcode
        try
        barBitmask = barcodeGen{idx}.rawBitmask;
        catch
        barBitmask = barcodeGen{idx}.bitmask;
        end
        
%         if isequal(sets.comparisonMethod,'dtw')
%             sets.idx = strcat(theoryStruct.name,num2str(idx));
%             comparisonStructure{idx} = ucr_dtw_score(theoryStruct.filename, barTested, barBitmask, sets);
%             comparisonStructure{idx}.bestBarStretch = 1;
%             comparisonStructure{idx}.length = lenBarTested;
%         else

        % run the loop for the stretch factors
        for j=1:length(stretchFactors)
            % here interpolate both barcode and bitmask 
            barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
            barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
            try
                [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or,rezMax{j}.secondPos,rezMax{j}.lengthMatch,rezMax{j}.dist] = comparisonFun(barC, theorBar, barB,theorBit,w);
            catch
                rezMax{j}.maxcoef = 0;rezMax{j}.pos=0;rezMax{j}.or=0;rezMax{j}.secondPos=0;rezMax{j}.lengthMatch=0;rezMax{j}.dist=0;
            end
            xcorrMax(j) = rezMax{j}.maxcoef(1);
        end


        % find which stretching parameter had the best score
        [value,b] = max(xcorrMax);


        % select the results for this best stretching parameter and output
        % them. If there were no values computed for this barcode, we don't
        % save anything.
        if ~isnan(value)
            rezMaxM{idx} = rezMax{b};
            bestBarStretch(idx) = stretchFactors(b);
            bestLength(idx) = round(lenBarTested*stretchFactors(b));
        end
%                 comparisonStructure{idx} = rezMax{b};
%                 comparisonStructure{idx}.bestBarStretch = stretchFactors(b);
%                 comparisonStructure{idx}.length = round(lenBarTested*stretchFactors(b));
%             else
%                  comparisonStructure{idx}.maxcoef(1:3) = nan;
%                  comparisonStructure{idx}.bestBarStretch = nan;
%                  comparisonStructure{idx}.length = nan;
%                  comparisonStructure{idx}.pos(1:3) = nan;
%                  comparisonStructure{idx}.or(1:3) = nan;
%             end
          % can also add ucr score here for convenience
%         if sets.comparison.useDTW
  %             % check that the positions are returned correctly, and how can
%             % this be implemented as an alternative to pcc, and how
%             % Sakoe-Chiba band corresponds to stretching
%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets); 
%         end
    end 

end
% % 
% %% test pccs
% pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
% bar1 = barTested; % or stretched if stretch factor is not 1
% bar2 = theorBar;
% bit1 = barBitmask;
% idxpos = rezMaxM{1}.secondPos(1);
% pos = rezMaxM{1}.pos(1);
% 
% if ~rezMaxM{1}.or
%     frag1 = flipud(bar1);
% else
%     frag1 = bar1;
% end
% % 
% % idxpos = idxpos+1
% frag2 = bar2(pos:pos+length(bar1)-1);
% % 
% bit1 = find(bit1,1,'first');
% % 
% idxpos = idxpos+bit1-1;
% frag1subseq = frag1(idxpos:idxpos+w-1);
% frag2subseq = frag2(idxpos:idxpos+w-1);
% % 
% pcc(frag1subseq',frag2subseq')
