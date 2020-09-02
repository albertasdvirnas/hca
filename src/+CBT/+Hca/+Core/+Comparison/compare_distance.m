function [rezMax,bestBarStretch,bestLength,rezMaxAll] = compare_distance(barcodeGen,theoryStruct, sets, consensusStructs )
    % compare_distance
    % Compares experiments to theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theoryStruct: theory structure
    %         consensusStructs: consensus structure
    % 
    %     Returns:
    %         comparisonStruct: comparison structure
        
    if nargin < 4
        consensusStructs = [];
    end
    disp('Starting comparing exp to theory...')
    tic
    % if there is consensus, add it as last element to barcodeGen
    if sets.genConsensus && ~isempty(consensusStructs)
        try
            for idx=1:length(consensusStructs) % add all possible consensuses
                barcodeGen{end+1} = consensusStructs{idx};
            end
        catch
            barcodeGen{end+1} = consensusStructs;
        end
    end
    
    stretchFactors = sets.theory.stretchFactors;
    comparisonMethod = sets.comparisonMethod;
%     islinear = sets.islinear
    w = sets.w;
    numPixelsAroundBestTheoryMask = 20;
            
    rezMax = cell(1,length(theoryStruct));
    bestBarStretch = cell(1,length(theoryStruct));
    bestLength = cell(1,length(theoryStruct));
    rezMaxAll = cell(1,length(theoryStruct));
%     rezMaxAll = 

    % Computing distances for each theory against all experiments. This
    % loop can be parallelized (parfor)
    parfor barNr = 1:length(theoryStruct)
        disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));
      
%         import CBT.Hca.Core.Comparison.on_compare_mp_all;
%         [rezMax{barNr},bestBarStretch{barNr},bestLength{barNr},rezMaxAll{barNr}] = on_compare_mp_all(barcodeGen,theoryStruct{barNr},comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask);

        import CBT.Hca.Core.Comparison.on_compare;
        [rezMax{barNr},bestBarStretch{barNr},bestLength{barNr}] = on_compare(barcodeGen,theoryStruct{barNr},comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask);
% % 
%         import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
%         comparisonStruct{barNr} = on_compare_theory_to_exp(barcodeGen,theoryStruct{barNr}, sets);
    end
    
      
    
    
    
    timePassed = toc;
    disp(strcat(['Experiments were compared to theory in ' num2str(timePassed) ' seconds']));
end

