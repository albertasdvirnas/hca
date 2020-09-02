function [dist] = unmasked_MASS_PCC(longVec, shortVec, shortVecBit,k)
    % masked MASS_PCC where short vector can have bitmask
    %
    %   Args:
    %       shortVec, longVec, shortVecBit,k
    %
    %   Returns:
    %
    %       dist
    
    shortVecCut = shortVec(logical(shortVecBit));

    % regular MASS PCC
    [dist] = MASS_PCC([longVec,longVec(1:sum(shortVecBit)-1)],shortVecCut, k);
    
    % shift back
    dist(1,:) = circshift(dist(1,:),[0,1-find(shortVecBit,1,'first')]);
    dist(2,:) = circshift(dist(2,:),[0,1-find(shortVecBit,1,'first')]);
    
    
end

