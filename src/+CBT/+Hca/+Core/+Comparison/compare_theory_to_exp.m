function [comparisonStruct] = compare_theory_to_exp(barcodeGen,theoryStruct, sets, consensusStructs )
    % compare_theory_to_exp
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
    if sets.genConsensus && ~isempty(consensusStructs)
        try
            for idx=1:length(consensusStructs) % add all possibleconsensuses
                barcodeGen{end+1} = consensusStructs{idx};
            end
        catch
            barcodeGen{end+1} = consensusStructs;
        end
    end
    
    stretchFactors = sets.theory.stretchFactors;
    comparisonMethod = sets.comparisonMethod;
    w = 200;
            
    rezMax = cell(1,length(theoryStruct));
    bestBarStretch = cell(1,length(theoryStruct));
    bestLength = cell(1,length(theoryStruct));

    % unfiltered comparison
    parfor barNr = 1:length(theoryStruct)
        disp(strcat(['comparing to theory barcode ' num2str(barNr) '_' theoryStruct{barNr}.filename] ));
        
%         import CBT.Hca.Core.Comparison.on_compare;
%         [rezMax{barNr},bestBarStretch{barNr},bestLength{barNr}] = on_compare(barcodeGen,theoryStruct{barNr},comparisonMethod,stretchFactors,w);
% 
        import CBT.Hca.Core.Comparison.on_compare_theory_to_exp;
        comparisonStruct{barNr} = on_compare_theory_to_exp(barcodeGen,theoryStruct{barNr}, sets);
    end
    
    timePassed = toc;
    disp(strcat(['Experiments were compared to theory in ' num2str(timePassed) ' seconds']));
end

