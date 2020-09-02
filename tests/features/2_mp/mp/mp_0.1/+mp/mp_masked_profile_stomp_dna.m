function [mp,mpI] = mp_masked_profile_stomp_dna(aTS, bTS, bitTS, r, kk)
    % mp for the case when bar1 is masked
    
    
    shortVecCut = aTS(logical(bitTS));

    
    import mp.mp_profile_stomp_dna;
    [mp,mpI] =  mp_profile_stomp_dna(shortVecCut, [bTS; bTS(1:r-1)],r,kk);
    
    % those who
    indexes = mpI<=length(bTS);
    % we substract back. We also need to substract the position on mpI!
    mpI(indexes) = mod(mpI(indexes)-find(bitTS,1,'first')-find(indexes)+1,length(bTS)) +1;
    
%     mpIF = flipud(mpI);

    indexes2 = mpI>length(bTS);
    
    % now we also transpose the indexes for the flipped version of barcode

    indx = [length(mpI):-1:1]';
     mpI(indexes2) = mod(mpI(indexes2) - indx(indexes2) -find(bitTS,1,'first')+1,length(bTS))+1;
     
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

