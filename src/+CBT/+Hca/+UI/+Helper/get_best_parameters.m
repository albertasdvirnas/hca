function [ maxcoef,pos,or ] = get_best_parameters( xcorrs, numBestPar, lenQ, isLinearTF, mask,theoryBit)
    % get_best_parameters
    
    % input xcorrs - cross correlation coefficients
    % numBestPar - number of coefficients to output
    % output rezMax
    
    if nargin < 2
        numBestPar = 3; % default is 3
    end
    
    if nargin < 3
       lenQ = 1;
    end
    
	if nargin < 4
       isLinearTF = 0;
    end
    
    if nargin < 5
        mask = 0; % default to previous if no mask
    end
    
    % invalid coefficients are treated as nan's.     
    if isLinearTF==1
        xcorrs(:,end-lenQ+2:1:end) = nan;
    end
    
    if nargin >=6
        % we can compute theoryBit for all positions along xcorrs
        nanValues = movmean(theoryBit,lenQ,'Endpoints','discard');
        xcorrs(:,nanValues > 50) = nan;        % 50 is estimate.., mean approx 50 bp  (1/10th) in a pixel were nan's
    end
    
    % now, based on lomgVecBit, bitmask some positions, maybe 
    %     longVecBit = [longVecBit  longVecBit(1:min(end-1,sum(shortVecBit)-1))];
    % shortVecCut

    % todo: based on options, don't record just the three best ones, but
    % introduce a "mask" for the best result, so that the same place is not reported
    % twice !

    % f gives the max over the two rows, s gives the indice of
    % which row has the max (1 or 2). TODO: maybe change this to be a match
    % table instead - corresponds to structural variations better.
    maxcoef = zeros(1,numBestPar);
    or = zeros(1,numBestPar);
    pos = zeros(1,numBestPar);

    for ii=1:numBestPar
        [f,s] =nanmax(xcorrs);

        % sort the max scores, ix stores the original indices
        [ b, ix ] = sort( f(:), 'descend','MissingPlacement','last' );

        % choose the best score
        indx = b(1) ;
        % save the best max score and orientation
        maxcoef(ii) = indx;
        or(ii) = s(ix(1));

        % finally, save the position. This can have two cases,
        % depending on the value of s
        pos(ii) = ix(1);
        
        % now add nan's to the area around the position of maxcoeff  
        xcorrs(:,max(1,ix(1)-mask):min(ix(1)+mask,size(xcorrs,2))) = nan;
    end
    
%     
%     
%            [f,s] =max(xcorrs);
% 
%         % sort the max scores, ix stores the original indices
%         [ b, ix ] = sort( f(:), 'descend' );
% 
%         % choose the best three scores. (change this in the future?)
%         indx = b(1:numBestPar)' ;
% 
%         % save the best three max score, and their orientation
%         maxcoef = indx;
%         or = s(ix(1:numBestPar)');
% 
%         % finally, save the position. This can have two cases,
%         % depending on the value of s
%         pos = ix(1:numBestPar)';

        
end
