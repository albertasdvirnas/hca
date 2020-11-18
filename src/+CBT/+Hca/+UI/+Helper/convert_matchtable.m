function [matchTable] = convert_matchtable(randStruct)
% convert data to query or query to tada match table


try
    matchtableOld = randStruct.matchTable;
catch
    matchtableOld = randStruct;

end
    matchTable = matchtableOld;


for i=1:size(matchtableOld,1)
    if matchtableOld(i,5) == 2;
            matchTable(i,1:4) = [matchtableOld(i,4)  matchtableOld(i,3) matchtableOld(i,2) matchtableOld(i,1)];
    else
        matchTable(i,1:4) = [ matchTable(i,3)  matchTable(i,4)  matchTable(i,1)  matchTable(i,2) ];
    end
end   


end

