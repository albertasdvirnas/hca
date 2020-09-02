function [maxSpearman,posMax,orMax,secondPos,lenM,scoreVec] = unmasked_spearman_corr2(A,B,W1)
% unmasked_spearman_corr - computes spearmancorrelation coefficient between
% A and B with lag and returns position and value of the best correlation
%

%   Args:
%       A - query
%       B - data
%       W1 - mask A (of type 000001111111100000000
%
%   Returns:
%
%     [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or,rezMax{j}.secondPos, rezMax{j}.abMatchTemp]

    Aw = compute_ranks(A(logical(W1)));
    % now convert Aw to rank
%     AwRank = compute_ranks(Aw);

    Awrev = flipud(Aw);
    maxSpearman = -inf;
    posMax = nan;
    orMax = 1;
    
    lenM = length(Aw);
    
    scoreVec = zeros(2,length(B)-length(Aw)+1);
    
    den = length(Aw)*(length(Aw).^2-1);
    for i=1:length(B)-length(Aw)+1
        Brank = compute_ranks(B((i:i+length(Aw)-1)));
        
        % could make this recursive if we fix ranks of B
        rhoF = 1-(6*sum((Aw-Brank).^2))/den;
        rhoRev = 1-(6*sum((Awrev-Brank).^2))/den;

%         [rhoF,~] = corr(Aw,B((i:i+length(Aw)-1)),'Type','Spearman');
%         [rhoRev,~] = corr(Awrev,B((i:i+length(Aw)-1)),'Type','Spearman');
        [score,orM] = max([ rhoF,rhoRev] );
        scoreVec(:,i) =  [rhoF;rhoRev];
        if maxSpearman < score
            maxSpearman = score;
            posMax = i;
            orMax = orM;
        end
            
    end

    % secondPos depends on the bitmask W1, and shows where on the A barcode
    % does the match start. Compare this to MP where it shows where
    % subfragment starts.
    secondPos = find(W1,1,'first');
    
    posMax = posMax-secondPos+1;
    


end

