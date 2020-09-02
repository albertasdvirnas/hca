function [ ] = export_pval_struct( rS,vals,data )

    fid = fopen(rS,'w');
    for i = 1:length(data)
        fprintf(fid, '%5d ', vals(i));
        fprintf(fid, '%5.4f ', data{i});
        fprintf(fid, '\n');
    end
    fclose(fid);


end

