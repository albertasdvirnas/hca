function [rezMax] = ucr_dtw_score(theory, shortVec, shortVecBit, sets)
    %   ucr_dtw_score - computes dtw score based on "trillion" code from UCR
    %
    %
    %   Args:
    %       shortVec, longVec, shortVecBit
    %   returns:
    %       rezMax - which stores maxcoef,pos, and or
    %
    
    shortVecCut = shortVec(logical(shortVecBit));

    % rand number, later change this to idx of barcode 
    nameFiles = sets.idx;
    % length of experiment
    M = sum(shortVecBit);
    % Sakoe-Chiba band, this corresponds to stretch factor
    R = sets.theory.stretchFactors(end)-1;

    % save experiment in temporary txt file. Check how this behaves in case
    % parfor is used
    fname1 = strcat([nameFiles 'query.txt']); fileID = fopen(fname1,'w');
    fprintf(fileID,'%2.5f ',shortVecCut); fclose(fileID);
    
    
    fname2 = strcat([nameFiles 'queryrev.txt']); fileID = fopen(fname2,'w');
    fprintf(fileID,'%2.5f ',fliplr(shortVecCut)); fclose(fileID);
    
    
    pathToScript = fullfile(pwd,sets.dtwscriptpath,'ucr_dtw.sh');

    outFile = fullfile(sets.output.matDirpath,strcat([nameFiles 'output.txt']));
%     tic
    ucrCode = fullfile(pwd,sets.dtwscriptpath,'a.out');

    cmdStr       = [pathToScript ' ' fname1 ' ' fname2 ' ' theory ' ' num2str(M) ' ' num2str(R) ' ' sets.output.matDirpath ' ' outFile ' ' ucrCode];
    system(cmdStr);
%     toc
    delete(fname1);
    delete(fname2);

    % smart would be to use mpi to save to different parts of the file..
    A = importdata(outFile);
    delete(outFile);
    
    % just make the coeff negative, so we look for max instead of min
    coef = -[A(2) A(4)];
    pos = [A(1) A(3)];
    [rezMax.maxcoef,rezMax.or ] = max(coef);
    rezMax.pos = pos(rezMax.or)+1-find(shortVecBit,1,'first')+1;
    
    
%     % mex the c function. This Should be mexed before    
%     mex 'bin/OVERLAPPING_DTW_MEX.cpp';
    
	% Instead of saving these every time, just have them saved as txt's and
	% pass them to ucr function. Should be a more clever way to do this
	% though.

%     % reversed
%     fname = strcat(['experimentrev.txt']); fileID = fopen(fname,'w');
%     fprintf(fileID,'%2.5f ',fliplr(shortVec(shortVecBit))); fclose(fileID);
%     
%     
%     % subsequence length.
%     N = sum(shortVecBit);
% 
%         % these first two are just to see how the query places as a
%     % subsequence, they're not comparable to the rest of the scores because
%     % they are for different length.
%     [pos(1,1), scores(1,1)] = OVERLAPPING_DTW_MEX(theory,'experiment.txt', N, sets.comparison.R);
%     [pos(1,2), scores(1,2)] = OVERLAPPING_DTW_MEX(theory,'experimentrev.txt', N, sets.comparison.R);

%     % then we have one of the two choices for position
%     if scores(1)< scores(2)
%         rezMax.maxcoef = scores(1);
%         rezMax.pos = pos(1);
%         rezMax.or = 1;
%     else
%         rezMax.maxcoef = scores(2);
%         rezMax.pos = pos(2);
%         rezMax.or = 2;   
%     end
    % faster function, only when barC has bitmask only on left and right
    %     xcorrs = unmasked_pcc_corr(barC, theorBar, barB);
    %     [rezMax.maxcoef,rezMax.pos,rezMax.or] = get_best_parameters(xcorrs, 3 );
    %     % now find the maximum score for this stretching parameter
    %     xcorrMax(j) = rezMax{j}.maxcoef(1);
                    
end

