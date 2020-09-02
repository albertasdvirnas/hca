% we run C script using bash function

sets.output.matDirpath = pwd;
query = imgaussfilt(rand(1,500),3);
data = imgaussfilt(rand(1,10000),3);


% query = 1:10;
% data = 1:10;
% data(4) = 5;

fname1 = strcat(['query.txt']); fileID = fopen(fname1,'w');
fprintf(fileID,'%2.5f ',query); fclose(fileID);
% reversed
fname2 = strcat(['queryrev.txt']); fileID = fopen(fname2,'w');
fprintf(fileID,'%2.5f ',fliplr(query)); fclose(fileID);



theory = strcat(['data.txt']); fileID = fopen(theory,'w');
fprintf(fileID,'%2.5f ',data); fclose(fileID);

M = length(query);
R = 0.01;

randNr=1;
% zscore(query)
sets.dtwscriptpath = '/home/albyback/git/hca/tests/features/6_dtw/Tests';
sets.dtwscriptpath = '/home/albyback/git/Development/hca/tests/features/6_dtw/Tests';

pathToScript = fullfile(sets.dtwscriptpath,'ucr_dtw.sh');
ucrCode = fullfile(sets.dtwscriptpath,'a.out');

outFile = fullfile(sets.output.matDirpath,strcat([num2str(randNr) 'output.txt']));
tic
cmdStr       = [pathToScript ' ' fname1 ' ' fname2 ' ' theory ' ' num2str(M) ' ' num2str(R) ' ' sets.output.matDirpath ' ' outFile ' ' ucrCode];
system(cmdStr);
toc
%% now test the result using matlab's in-built DTW at best position

A= importdata(outFile);
A(2)
A(1)

MAXSAMP = round(M*R);
pos = A(1)+1;
bestPosData = data(pos:pos+M-1);
% [DIST,IX,IY] = dtw(zscore(bestPosData,1),zscore(query,1),MAXSAMP);
% [DIST,IX,IY] = dtw(zscore(query,1),zscore(bestPosData,1),'squared',MAXSAMP-2);
[DIST,IX,IY] = dtw(zscore(query,1),zscore(bestPosData,1),'squared',MAXSAMP);
% DIST
sqrt(DIST)

% sqrt(sum((zscore(query,1)-zscore(bestPosData,1)).^2))
