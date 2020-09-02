function [transformedValues] = compute_ranks(B)

    % transformedValues
    % convert values of A to ranks
    values = 1:length(B);
    % sorted values
    [~,b] = sort(B);
    
    transformedValues = zeros(length(B),1);

    transformedValues(b) = values;

end

