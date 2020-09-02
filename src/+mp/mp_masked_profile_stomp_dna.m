function [mp,mpI,indexes] = mp_masked_profile_stomp_dna(aTS, bTS, bitTS, r, kk,islinear)
    % mp for the case when bar1 is masked
    
    % Args:
    
    % Returns:
    %
    % mp - matrix profile
    % mpI - matrix profile index
    % indexes - orientation of matrix profile indexes (forward 1 reverse 0)
    
    if nargin < 6
        islinear = 0;
    end
    
    shortVecCut = aTS(logical(bitTS));

    %     numPos - number of positions for forward/reverse match
    import mp.mp_profile_stomp_dna;
    if islinear
        [mp,mpI] =  mp_profile_stomp_dna(shortVecCut, bTS,r,kk);
        numPos = length(bTS)-r+1;
    else
        [mp,mpI] =  mp_profile_stomp_dna(shortVecCut, [bTS; bTS(1:r-1)],r,kk);
        numPos = length(bTS);
    end
    
    % forward indices,
    indexes = mpI<=numPos;
    % we substract back. We also need to substract the position on mpI! But
    % mpI still keeps information about this position
%     mpI(indexes) = mod(mpI(indexes)-find(bitTS,1,'first')-find(indexes)+1,length(bTS)) +1;
     mpI(indexes) = mpI(indexes)-find(bitTS,1,'first')-find(indexes)+2;
    
%     mpIF = flipud(mpI);

    % backward indices
    indexes2 = mpI>numPos;
    
    % now we also transpose the indexes for the flipped version of barcode,
    % because for the reverse these are inverted(?) are they 
    indx = [length(mpI):-1:1]';
    
%      mpI(indexes2) = mod(mpI(indexes2) - indx(indexes2) -find(bitTS,1,'first')+1,length(bTS))+1;
   mpI(indexes2) = mpI(indexes2) - indx(indexes2) -find(bitTS,1,'first')+2-numPos;
     
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

