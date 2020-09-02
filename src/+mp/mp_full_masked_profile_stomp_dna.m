function [mp,mpI,indexes] = mp_full_masked_profile_stomp_dna(aTS, bTS, aBitTS,bBitTS, r, kk,islinear)
    % mp for the case when bar1 is masked and bar2 is masked
    
    % Args:
    
    % Returns:
    %
    % mp - matrix profile
    % mpI - matrix profile index
    % indexes - orientation of matrix profile indexes (forward 1 reverse 0)
    
    if nargin < 6
        islinear = 0;
    end
    
    % i think this should always be linear, since both of these have
    % bitmasks, no point in looping around to find the best match (could do
    % that in post processing if it was relevant..)
    
    shortVecCut = aTS(logical(aBitTS));
    % mask the bTS too
    bTS = bTS(logical(bBitTS));

    %     numPos - number of positions for forward/reverse match
    import mp.mp_profile_stomp_dna;
    if islinear
        [mp,mpI] =  mp_profile_stomp_dna(shortVecCut, bTS,r,kk);
        numPos = length(bTS)-r+1;
    else
        [mp,mpI] =  mp_profile_stomp_dna(shortVecCut, [bTS; bTS(1:r-1)],r,kk);
        numPos = length(bTS);
    end
    
    %%
%      [a,b] = max(mp);
%     b2 = mpI(b);
%     A = shortVecCut(b:b+r-1);
%     if b2 > numPos
%         B = bTS(b2-numPos:b2-numPos+r-1);
%     else
%         B = bTS(b2:b2+r-1);
%     end
%     pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
%     corr = pcc(A,B)
%     
%     % then we want the position on B converted to the start on A.
%     % this seems to be ok?
%     b2-numPos-b+1
%     
%     stA = b;
%     stB = b2-numPos;
%     
%     stA1 =  b2-numPos-b+1
    

    %%
    % forward indices,
    indexes = mpI<=numPos;
    indexes2 = mpI>numPos;

    % we substract back. We also need to substract the position on mpI! But
    % mpI still keeps information about this position
%     mpI(indexes) = mod(mpI(indexes)-find(bitTS,1,'first')-find(indexes)+1,length(bTS)) +1;
     mpI(indexes) = mpI(indexes)-find(aBitTS,1,'first')+find(bBitTS,1,'first')-find(indexes)+1;
%     find(bBitTS,1,'first')
%     mpIF = flipud(mpI);

    % backward indices
    
    % now we also transpose the indexes for the flipped version of barcode,
    % because for the reverse these are inverted(?) are they 
%     indx = [length(mpI):-1:1]';
     indx = [1:length(mpI)]';
    
%      mpI(indexes2) = mod(mpI(indexes2) - indx(indexes2) -find(bitTS,1,'first')+1,length(bTS))+1;
    mpI(indexes2) = mpI(indexes2) - indx(indexes2) -find(aBitTS,1,'first')+find(bBitTS,1,'first')+1-numPos;
%      mpI(262)-numPos-indx(262)+3-find(aBitTS,1,'first')-find(bBitTS,1,'first')
   
   %%
%      [a,b] = max(mp);
%     b2 = mpI(b);
%     A = shortVecCut(b:b+r-1);
% %     if b2 > numPos
%     B = bTS(b2+b+numPos:b2-numPos+r-1);
% %     end
%     pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
%     corr = pcc(flipud(A),B)
%     mpI(147)-length(bTS)-indx(147)-find(bitTS,1,'first')+2
% 
%     % mpI(1)-length(bTS)
%     % but this is 
%     pos =  mod( 2*length(bTS) - mpI(indexes2)-find(indexes2)+find(bitTS,1,'first')-3,length(bTS))+1;
% %     figure,plot()
% %   mpI(indexes2) =  length(bTS) - mpI(indexes2)
%   mod( 2*length(bTS) - mpI(indexes2)-find(indexes2)+find(bitTS,1,'first')-3,length(bTS))+1;
% %   
%   test = mod( 2*length(bTS)-mpI(indexes2)+find(bitTS,1,'first')+3-find(indexes2),length(bTS))+1;
%   test
%   2*length(bTS) - mpI(indexes2)+find(bitTS,1,'first')-3-find(indexes2)
%   
%   2*length(bTS) - mpI(indexes2)-find(indexes2)
% %      2*length(bTS)-mpI(1)-find(bitTS,1,'first')+1+1+1
%       2*length(bTS)-mpI(2)-find(bitTS,1,'first')+1+1
%             2*length(bTS)-mpI(3)-find(bitTS,1,'first')+1+1-1

%     mpd=
%     end/2
%     mpI(1:end/2) = mod(mpI(1:end/2)-find(bitTS,1,'first')+1);
%         mpI(end/2) = mpI(1:end/2)-find(bitTS,1,'first')+1;
% 
%     % shift back
%     dist(1,:) = circshift(dist(1,:),[0,1-find(shortVecBit,1,'first')]);
%     dist(2,:) = circshift(dist(2,:),[0,1-find(shortVecBit,1,'first')]);
%     
end

