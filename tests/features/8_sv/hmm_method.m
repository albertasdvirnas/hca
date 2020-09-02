% generate random barcode
len1 = 500;
len2 = 1000;   
sigma = 2.3;
pcc = 0.7;
name = 'test1';
[bar1,bar2,randPos] = generate_rand(len1, len2, pcc, sigma, islinear);

barStruct.bar1 = bar1;
barStruct.bar2 = bar2;


sets = ini2struct('sv_settings.txt');
if ~exist(sets.fold,'dir' )
    mkdir(sets.fold );
end
for i=1:length(barStruct)
    folder = strcat([sets.fold sets.svType '_' num2str(length(barStruct.bar1 )) '_' name '_single_possible']);

    % save DATA
    rS = strcat([folder num2str(i) '_data_seq.txt']);
    fid = fopen(rS,'w');
    fprintf(fid, '%5.4f ', barStruct.bar2  );
    fprintf(fid, '\n');
    fclose(fid);
    barStruct.name1 = rS;


    % save QUERY (this includes the variation)
    rS = strcat([folder num2str(i) '_query_seq.txt']);
    fid = fopen(rS,'w');
    fprintf(fid, '%5.4f ',barStruct.bar1  );
    fprintf(fid, '\n');
    fclose(fid);
    barStruct.name2 = rS;
end
%

sets.cg = 0;

import Core.run_hmm_matlab;
output2 = run_hmm_matlab({barStruct},sets,'test1');

tic
import Core.run_hmm_bash;
output = run_hmm_bash({barStruct},sets,name);
toc


