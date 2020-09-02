function [ combinedStruct ] = combine_chromosome_results(theoryStruct, comparisonStruct )
% combine_chromosome_results
%
% This function stores the top placement for each experiment
%
%     Args:
%         theoryStruct, comparisonStruct
% 
%     Returns:
%         combinedStruct: Return structure

    indx = zeros(1,length(rezMax{1}));
    combinedStruct = cell(1, length(comparisonStruct{1}));

    for i=1:length(comparisonStruct{1})
        try
            maxCoefs = cellfun(@(x) x{i}.maxcoef(1),comparisonStruct);
            [~,indx(i)] = max(maxCoefs);
            combinedStruct{i} = comparisonStruct{indx(i)}{i};
            combinedStruct{i}.idx = indx(i);
            combinedStruct{i}.pos =  comparisonStruct{indx(i)}{i}.pos; %+sum(accuLengths(1:indx(i)-1));
            combinedStruct{i}.name = theoryStruct{indx(i)}.name;
        catch
            disp(strcat(['No placement found for barcode ' num2str(i)]));
        end
        
    end

end

