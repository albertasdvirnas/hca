function [maxSpearman,posMax,orMax,secondPos,lenM,scoreVec] = unmasked_spearman_corr3(A,B,W1)
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
    
    [a,b] = sort(B((1:1+length(Aw)-1)));

    Brank = compute_ranks(B((1:1+length(Aw)-1)));
%     tic
%     BrankTotal = compute_ranks(B);
%     toc
    den = length(Aw)*(length(Aw).^2-1);
    
    % could make this recursive if we fix ranks of B
    rhoF = 1-(6*sum((Aw-Brank).^2))/den;
    rhoRev = 1-(6*sum((Awrev-Brank).^2))/den;

    [maxSpearman,orMax] = max([ rhoF,rhoRev] );
    posMax = 1;

%     values = 1:length(Brank);
    % sorted values
%     [~,b] = sort(B);
    
%     transformedValues = zeros(length(B),1);

%     transformedValues(b) = values;

    
        
    for i=2:length(B)-length(Aw)+1
        % now we remove the first element from Brank and add new element,
        % so update Brank
%         [~,b] = sort(B(i:i+length(Aw)-1));
%         Brank(b) = values;
        c = find(B(i+length(Aw)-1)<a,1,'first');
        Brank(b(1):end) = Brank(b(1):end) - 1;
        Brank = circshift(Brank,[1,0]);
        if ~isempty(c)
            Brank(end) = c;
            Brank(b(c):end) = Brank(b(c):end)+1;
        else
             Brank(end) = length(B);
        end
%         Brank = compute_ranks(B(i:i+length(Aw)-1));
        
%         Brank = compute_ranks(B((i:i+length(Aw)-1)));
%         Brank = circshift(Brank,[1,0])
%         valsMore = Brank > Brank(end);
%         Brank(valsMore) = Brank(valsMore)-1;
%         %%
%         Brank(1) = B(i+length(Aw)-1)

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

