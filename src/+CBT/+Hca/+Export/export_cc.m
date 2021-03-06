function [matFilepath] = export_cc(T, matDirpath) 
    % Exports table as txt

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([ 'table_' timestamp  '.txt']);
    matFilename2 = strcat([ 'table_' timestamp  '.dat']);


    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = strcat([matDirpath, matFilename]);
    writetable(T,matFilepath,'Delimiter','\t')  
    matFilepath = strcat([matDirpath, matFilename2]);

    writetable(T,matFilepath,'WriteRowNames',true)  
    fprintf('Saved table to ''%s''\n', matFilepath);
end
