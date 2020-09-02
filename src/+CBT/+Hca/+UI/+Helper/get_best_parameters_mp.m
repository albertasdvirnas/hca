function [ maxcoef,pos,or,idxpos ] = get_best_parameters_mp( mp,mpI, mpD, numBestPar, mask,islinear,lenD)
    % get_best_parameters_mp
    %
    % todo: merge to get_best_parameters, so that all cases are in the same
    %
    %   Args:
    %
    %   Returns:
    %
    % input xcorrs - cross correlation coefficients
    % numBestPar - number of coefficients to output
    % output rezMax
    
    if nargin < 4
        numBestPar = 3; % default is 3
    end
    
%     if nargin < 5
%        lenQ = 1;
%     end
%     
% 	if nargin < 6
%        isLinearTF = 0;
%     end
    
    if nargin < 5
        mask = 0; % default to previous if no mask
    end
    
    % there is no invalid coefficients anymore! in this newer version!
%     % invalid coefficients are treated as nan's.     
%     if isLinearTF==1
%         % valid indices 1..lenD-r+1
%         mpIvals = mpI >=lenD-r+2;
%         mpivalsD = mpI <=lenD;
%         mp(logical(mpIvals.*mpivalsD)) = nan;
% 
%         mpIvals = mpI >=2*lenD-r+2;
%         mpI(mpIvals) = nan;
% 
% %         mpIvals = mpI >=lenD-lenQ+2;
% %         mpivalsD = mpI <lenD;
% %         mp(logical(mpIvals.*mpivalsD)) = nan;
% %            
% %         mpIvals = mpI >=2*lenD-lenQ+2;
% %         mpI(mpIvals) = nan;
%      
% %         % find which mpI values are going over..
% %         xcorrs(:,end-lenQ+2:1:end) = nan;
%     end
    % todo: based on options, don't record just the three best ones, but
    % introduce a "mask" for the best result, so that the same place is not reported
    % twice !

    % f gives the max over the two rows, s gives the indice of
    % which row has the max (1 or 2). TODO: maybe change this to be a match
    % table instead - corresponds to structural variations better.
    maxcoef = zeros(1,numBestPar);
    or = zeros(1,numBestPar);
    pos = zeros(1,numBestPar);
    idxpos = zeros(1,numBestPar);

    for ii=1:numBestPar
        [f,s] =nanmax(mp);

        % sort the max scores, ix stores the original indices
        [ b, ix ] = sort( f(:), 'descend','MissingPlacement','last' );

        % choose the best score
        indx = b(1) ;
        % save the best max score and orientation
        maxcoef(ii) = indx;
        
        if islinear
            pos(ii) = mpI(s(ix(1)));
        else
            pos(ii) = mod(mpI(s(ix(1)))-1,lenD)+1;
        end
        
        % now orientation is easy!
        or(ii) = mpD(s(ix(1)));
        %double(mpI(s(ix(1)))>lenD)+1;
        
        idxpos(ii) = s(ix(1)); % should also add bitmask shift..

        % finally, save the position. This can have two cases,
        % depending on the value of s
        mp(s) = nan;
        
        % special case is when we are at the border of mpI, so that mask
        % would take to the other side in circular case..
        mp(logical((mpI > mpI(s(ix(1)))-mask).*(mpI < mpI(s(ix(1)))+mask))) = nan;
        
        % find values within these. take into account edges..
%         mp(logical((mod(mpI-1,lenD)+1 > pos(ii)-mask).*( mod(mpI-1,lenD)+1 < pos(ii)+mask))) = nan;
        % now have find values on mpI close to this
        % now add nan's to the area around the position of maxcoeff  
%         xcorrs(:,max(1,ix(1)-mask):min(ix(1)+mask,size(xcorrs,2))) = nan;
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