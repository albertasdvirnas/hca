function [ comparisonStructure ] = on_compare_theory_to_exp( barcodeGen,theoryStruct, sets)
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
    comparisonStructure = cell(1, length(barcodeGen));
    %ccM = cell(length(rawBarcodes),1);
    
    switch sets.comparisonMethod
        case 'unmasked_pcc_corr'
            import SignalRegistration.unmasked_pcc_corr;
            comparisonFun = @(x,y,z,w) unmasked_pcc_corr(x,y,z);
            import CBT.Hca.UI.Helper.get_best_parameters;

%             parameterfun = @(x)
        case 'mass_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));
        case 'dtw'
            
            import SignalRegistration.ucr_dtw_score;
            comparisonFun = @(x,y,z) ucr_dtw_score(x,y,z, sets);

%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets);
        case 'mp'
            import mp.mp_dist_stomp_with_masks;
            [scoreMatPCC,scorePos,orientation,secondPos, abMatchTemp] = mp_dist_stomp_with_masks(X1,X2, bitX1, bitX2, w, kk,par,islinear)

            comparisonFun = @(x,y,z,w) unmasked_MP(y,x,z,2^(4+nextpow2(length(x))));

            error('not yet implemented');
            
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
    stretchFactors = sets.theory.stretchFactors;
    try % some versions of sets might not have this number, therefore we add a test here
        minLength = sets.comparison.minLength;
    catch
        minLength = 20;
    end
    
    % number of pixels around best theory match coefficient
    try
        numPixelsAroundBestTheoryMask = sets.bitmasking.numPixels;
    catch
        numPixelsAroundBestTheoryMask = 0;
    end

    
    % load theory barcode txt file. For UCR DTW (c++ code), we only need the name of the
    % file so this can be skipped.
    fileID = fopen(theoryStruct.filename,'r');
    formatSpec = '%f';
    theorBar = transpose(fscanf(fileID,formatSpec));
    fclose(fileID);
%     theorBit = ones(1,length(theorBar));

   	import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode
% 	import SignalRegistration.unmasked_pcc_corr;

	%import SignalRegistration.masked_pcc_corr;
    
    % for all the barcodes run
    for idx=1:length(barcodeGen)

        % xcorrMax stores the  maximum coefficients
        xcorrMax = zeros(1,length(stretchFactors));
        
        % rezMaz stores the results for one barcode
        rezMax = cell(1,length(stretchFactors));
       
        % barTested barcode to be tested
        barTested = barcodeGen{idx}.rawBarcode;
        
        % in case barcode should be filtered
        barTested = filter_barcode(barTested, sets.filterSettings);

        % length of this barcode
        lenBarTested = length(barTested);
        
        % barBitmask - bitmask of this barcode
        barBitmask = barcodeGen{idx}.rawBitmask;
        
        if isequal(sets.comparisonMethod,'dtw')
            sets.idx = strcat(theoryStruct.name,num2str(idx));
            comparisonStructure{idx} = ucr_dtw_score(theoryStruct.filename, barTested, barBitmask, sets);
            comparisonStructure{idx}.bestBarStretch = 1;
            comparisonStructure{idx}.length = lenBarTested;
        else

            % run the loop for the stretch factors
            for j=1:length(stretchFactors)
                % TODO: choose interpolation method. I.e. could use the ideal
                % sinc interpolation for possibly better results

                % here interpolate both barcode and bitmask 
                barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
                barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));

                % compute the scores.
                % we limit here to computing theory for short experiment vs. 
                % long theory. In case the interpolated experiment is longer
                % than the theory, we display a warning, and output 0's for
                % this comparison
                % 
                if length(barC) > length(theorBar) % in case exp larger than theory
                    disp(strcat(['Warning, the experiment '  barcodeGen{idx}.name ' with stretching ' num2str(stretchFactors(j)) ' is longer than theory']))
                    xcorrMax(j) = nan;
                    %xcorrs = masked_pcc_corr(theorBar, barC, theorBit,barB);
                else
                    if sum(barB) <= minLength % move this to settings
                        disp(strcat(['Warning, after bitmasking, the experiment '  barcodeGen{idx}.name ' with stretching ' num2str(stretchFactors(j)) ' is shorter than 20 pixels']))
                        xcorrMax(j) = nan;    
                    else
                        % faster function, only when barC has bitmask only on left and right
                        xcorrs = comparisonFun(barC, theorBar, barB);
                        [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or] = get_best_parameters(xcorrs, 3,length(barC),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);
                        % now find the maximum score for this stretching parameter
                        xcorrMax(j) = rezMax{j}.maxcoef(1);
                    end
                end             
            end       


            % find which stretching parameter had the best score
            [value,b] = max(xcorrMax);



            % select the results for this best stretching parameter and output
            % them. If there were no values computed for this barcode, we don't
            % save anything.
            if ~isnan(value)
                comparisonStructure{idx} = rezMax{b};
                comparisonStructure{idx}.bestBarStretch = stretchFactors(b);
                comparisonStructure{idx}.length = round(lenBarTested*stretchFactors(b));
            else
                 comparisonStructure{idx}.maxcoef(1:3) = nan;
                 comparisonStructure{idx}.bestBarStretch = nan;
                 comparisonStructure{idx}.length = nan;
                 comparisonStructure{idx}.pos(1:3) = nan;
                 comparisonStructure{idx}.or(1:3) = nan;
            end
        end
          % can also add ucr score here for convenience
%         if sets.comparison.useDTW
  %             % check that the positions are returned correctly, and how can
%             % this be implemented as an alternative to pcc, and how
%             % Sakoe-Chiba band corresponds to stretching
%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets); 
%         end
    end 

end

