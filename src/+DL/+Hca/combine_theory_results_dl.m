function [ rezMaxC ] = combine_theory_results_dl(theoryStruct, rezMax)
% combine_theory_results_dl
%
% This function stores the top placement for each experiment and each
% barcode type
%
%     Args:
%         theoryStruct, comparisonStruct
%
%     Returns:
%         combinedStruct: Return structure

rezMaxC = cell(1, length(rezMax{1}));

for i=1:length(rezMax{1})
  for btype = ["dual", "dense", "sparse"]
    thisRezMax = cellfun(@(x) ...
      cellfun(@(y) subsref(y, substruct('.', btype)), x, 'un', 0), ...
      rezMax, 'un', 0);
    maxCoefs = cellfun(@(x) x{i}.maxcoef(1), thisRezMax);
    [~, thrInd] = nanmax(maxCoefs);
    newRexMax = thisRezMax{thrInd}{i};
    newRexMax.idx = thrInd;
    newRexMax.name = theoryStruct{thrInd}.name;
    rezMaxC{i} = subsasgn(rezMaxC{i}, substruct('.', btype), newRexMax);
  end
end
