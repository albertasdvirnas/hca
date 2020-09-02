function outStruct = save_rand(bar1,folder)
    %
    %
    %
    outStruct = [];
    for i=1:length(bar1)

   
         % save the variation in the folder
        rS = strcat([folder num2str(i) '_seq.txt']);
        fid = fopen(rS,'w');
        fprintf(fid, '%5.4f ',bar1{i} );
        fprintf(fid, '\n');
        fclose(fid);
        outStruct{i}.filename = rS;
%         rS = strcat([folder num2str(i) '_query_seq.txt']);
%         fid = fopen(rS,'w');
%         fprintf(fid, '%5.4f ', structBar{i}.bar2 );
%         fprintf(fid, '\n');
%         fclose(fid);
%         rS = strcat([folder num2str(i) '_pos_seq.txt']);
%         fid = fopen(rS,'w');
%         fprintf(fid, '%5.4f ', pos(1,:) );
%         fprintf(fid, '\n');
%         fprintf(fid, '%5.4f ', pos(2,:));
%         fprintf(fid, '\n');
%         fclose(fid);
%         
%         rS = strcat([folder num2str(i) '_matchtable.txt']);
%         fid = fopen(rS,'w');
%         for j=1:size(matchTable,1)
%             fprintf(fid, '%4d ', matchTable(j,:) );
%             fprintf(fid, '\n');
%         end
%         fclose(fid);
%         
        % add noise
%         import Rand.add_normal_noise;
%         import Rand.add_noise_to_barcode;
% 
%         [bar2] = add_noise_to_barcode(bar2, sets);

%         % can also save discretized series/fasta files
%         import Comparison.discretize_series;
%         [ bar1D ] = discretize_series( zscore(bar1), sets.alfabet_size);
%         [ bar2D ] = discretize_series( zscore(bar2), sets.alfabet_size);
% 
%         
%         % convert to protein
%         dataAA = int2aa(bar1D);
%         queryAA = int2aa(bar2D);

%         % convert to
%         rS = strcat([folder num2str(i) '_query_prot.fasta']);
%         fd = fopen(rS,'w');
%         fprintf(fd,strcat(['>'  'query_' num2str(i) '\n']));  
%         fprintf(fd,'%s', structBar{i}.queryAA);  
%         fprintf(fd,'\n'); 
%         fclose(fd);
%         rS = strcat([folder num2str(i) '_data_prot.fasta']);
%         fd = fopen(rS,'w');
%         fprintf(fd,strcat(['>'  'data_' num2str(i) '\n']));  
%         fprintf(fd,'%s', structBar{i}.dataAA);  
%         fprintf(fd,'\n'); 
%         fclose(fd);
    end
end

