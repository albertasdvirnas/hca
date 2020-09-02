function [corr] = mp_corr_test2(bar1, bar2,bit1, bit2, r,islinear)



% bar2 = [imgaussfilt(normrnd(0,1,len2,1),2.3); bar1; imgaussfilt(normrnd(0,1,len2,1),2.3)];

% bit1 = ones(length(bar1),1);

% r = 300;
k= 2^14;
% isLinearTF = 0;

% this should compute linear as well (because otherwise we compute too much
import mp.mp_dist_stomp_with_masks;
[maxcoef,pos,or,idxpos ] = mp_dist_stomp_with_masks(bar1, bar2,bit1, bit2, r, k,inf);
% 
% import mp.mp_masked_profile_stomp_dna;
% comparisonFun = @(x,y,z,w,u) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))),u);
% 
% tic
% [mp, mpI,mpD] = comparisonFun(bar1,bar2,bit1,r,islinear);
% toc
% 
% % what we want:
% % 1) to have placement of (whole) bar1 along bar2
% % 2) to be able to check that the placement of subsequence of bar1 along
% % bar2 has correct dist score as described in mp
% % Finally run through the functions to see
% 
% import CBT.Hca.UI.Helper.get_best_parameters_mp;
% [ maxcoef,pos,or,idxpos ] = get_best_parameters_mp( mp,mpI,mpD, 1, 50,islinear,length(bar2));
% % or
pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);

if ~or(1)
    frag1 = flipud(bar1);
else
    frag1 = bar1;
end

% frag2 = bar2(pos(1):pos(1)+length(bar1)-1);

bit1 = find(bit1,1,'first');

% ok for this example
% frag1(idxpos-bit1+1:idxpos-bit1+1+r-1);
frag1subseq = frag1(idxpos+bit1-1:idxpos+bit1-1+r-1);
% this should be from 1..
% frag2subseq = bar2(pos(1)+idxpos-bit1:pos(1)+idxpos-bit1+r-1);

frag2subseq = bar2(pos(1)+idxpos+bit1-2:pos(1)+idxpos+bit1+r-3);

corr = pcc(frag1subseq,frag2subseq);
% ok!



end

