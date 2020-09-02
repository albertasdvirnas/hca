function [scoreMatPCC,scorePos,orientation,secondPos, w,dist] = mp_dist_stomp_with_masks(X1,X2, bitX1, bitX2, w, kk,par,islinear)

    %   Args:
    %           X1
    %           X2
    %           w
    %
    %   Returns:
    %       mpd
    
    if nargin < 8
        islinear = 0;
    end
    
    islinear = 1;
    dist = [];
%     n = length(X1);
%     m = length(X2);
    
    % in this case both barcodes can be bitmasked, i.e. both X1 and X2. We
    % treat this as linear case, i.e. we do not circularly shift around
    % either of the barcodes (because if we did then bitmasked region would
    % appear in the middle of the barcode sometimes, which does not make
    % sense)
    
%     import mp.mp_profile_stomp_dna;
    
    import mp.mp_full_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z1,z2,w,u) mp_full_masked_profile_stomp_dna(x,y,z1,z2,w,kk,u);


    
    % use bitmasks, at the moment they are added here, but not used in
    % calculation (so are assumed all to be ones)
    %     X1 = X1(logical(bitX1));
    %     X2 = X2(logical(bitX2));
    
    % compute, non-circular
    % P_AB
    % X1 is allowed to have a bitmask
    [mpAB,mpIAB,mpDAB] = comparisonFun(X1, X2, bitX1, bitX2, w,islinear);
    
    % now if we want to plot the fragments and
    
    % allow circular in getting best parameters, but it might be negative
    import CBT.Hca.UI.Helper.get_best_parameters_mp;
    [ maxcoef,pos,or,idxpos ] = get_best_parameters_mp( mpAB,mpIAB,mpDAB, 1, 50, islinear, length(X2));
% or
    % TODO: check if islinear = 1 works with flipped orientation
% 
%         if ~or(1)
%             frag1 = flipud(X1);
%         else
%             frag1 = X1;
%         end
% 
%         frag2 = bar2(pos(1):pos(1)+length(bar1)-1);
%         % we get wrong value as output sometime?
%         bit1 = find(bitX1,1,'first');
%         bit2 = find(bitX2,1,'first');
% 
%         % ok for this example
%         % frag1(idxpos-bit1+1:idxpos-bit1+1+r-1);
%         frag1subseq = frag1(idxpos+bit1-1:idxpos+bit1-1+w-1);
%         % this should be from 1..
%         % frag2subseq = bar2(pos(1)+idxpos-bit1:pos(1)+idxpos-bit1+r-1);
% 
%         frag2subseq = X2(pos(1)+idxpos+bit1-2:pos(1)+idxpos+bit1+w-3);
% 
%         corr = pcc(frag1subseq,frag2subseq);
% ok!

    scoreMatPCC = maxcoef;
    scorePos = pos;
    orientation = or;
    secondPos = idxpos;
% %     
% %     %%%TEST
% % 
%     pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
% 
% 
%     % frag2 = bar2(pos(1):pos(1)+length(bar1)-1);
% 
%     bit1 = find(bitX1,1,'first');
%     bit2 = find(bitX2,1,'first');
% % % 
%     dd = 0;
%     % ok for this example
%     frag1subseq = X1(idxpos+bit1-1:idxpos+bit1-1+w-1);
%     
% % %     
%     if ~or(1)
%         frag1subseq = flipud(frag1subseq);
%     end
% % %     
% % %     % this should be from 1..
%     frag2subseq = X2(pos(1)+idxpos+bit1-1+dd:pos(1)+idxpos+bit1-1+w-1+dd);
% % % 
%     corr = pcc(frag1subseq,frag2subseq)


% ok!
    %     indexes = mpIAB <= length(X2)-w+1;
    %     % we substract back. We also need to substract the position on mpI!
    %     mpIAB(indexes) = mod(mpIAB(indexes)-find(bitX1,1,'first')-find(indexes)+1,length(X2)-w+1) +1;
    %     
    %     indexes2 = mpIAB>length(X2)-w+1;

    %     % now we also transpose the indexes for the flipped version of barcode
    %     indx = [length(mpIAB):-1:1]';
    %     mpIAB(indexes2) = mod(mpIAB(indexes2) - indx(indexes2) -find(bitX1,1,'first')+1,length(X2)-w+1)+1;
    %      
     
    %% IF we also do mpBA: but as result we should always outbut bar1 on bar2 comparison
%     abMatchTemp = 1;
% 
%     % P_BA
%     [mpBA,mpIBA,mpDBA] = comparisonFun(X2, X1,bitX2,bitX1, w, islinear);
%     
%        
%     % allow circular in getting best parameters, but it might be negative
%     import CBT.Hca.UI.Helper.get_best_parameters_mp;
%     [ maxcoef2,pos2,or2,idxpos2 ] = get_best_parameters_mp( mpBA,mpIBA,mpDBA, 1, 50, islinear, length(X1));
% or

%    indexes = mpIBA <= length(X1)-w+1;
%     % we substract back. We also need to substract the position on mpI!
%     mpIBA(indexes) = mod(mpIBA(indexes)-find(bitX2,1,'first')-find(indexes)+1,length(X1)-w+1) +1;    
%     indexes2 = mpIBA>length(X1)-w+1;   
%     % now we also transpose the indexes for the flipped version of barcode
%     indx = [length(mpIBA):-1:1]';
%     mpIBA(indexes2) = mod(mpIBA(indexes2) - indx(indexes2) -find(bitX2,1,'first')+1,length(X1)-w+1)+1;
   
%     % P_ABBA
%     mpABBA = [mpAB;mpBA];
%     % I_AMMA
%     mpIABBA = [mpIAB;mpIBA];
% 
%     % k, specific value.
%     % Don't take the highest value, but the one close to the top values..
%     k = par*(n+m);
% 
%     if length(mpABBA) > k
%         [a,~] = sort(mpABBA,'desc');
%         mpd = a(round(k));
%         [~,mpdI] = max(mpABBA);
%         mpdII = mpIABBA(mpdI);
%     else
%         [mpd,mpdI] = max(mpABBA); % max, gives the highest pcc as distance between the two, also would be nice to have position
%         mpdII = mpIABBA(mpdI);
%     end
%     
    %%
%     X1(mpdI:mpdI+w-1)





end

