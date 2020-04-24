function [outputArg1,outputArg2] = unmasked_dtw_score(theory, shortVec,shortVecBit)


    import SignalRegistration.ucr_dtw_score;
    [rezMax] = ucr_dtw_score(theory, experiment, expBit, sets)
    shortVecCut = shortVec(logical(shortVecBit));

    % regular MASS PCC
    [dist] = MASS_PCC([longVec,longVec(1:sum(shortVecBit)-1)],shortVecCut, k);
    
    % shift back
    dist(1,:) = circshift(dist(1,:),[0,1-find(shortVecBit,1,'first')]);
    dist(2,:) = circshift(dist(2,:),[0,1-find(shortVecBit,1,'first')]);
    
end

