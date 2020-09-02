function [pos] = mp_position_test(len2,bitmask,islinear)
    % run_mp_generation_tests
    % MP generation tests to see if Matrix Profile method runs properly (from both computational 
    % and results perspective ) for
    % sequence similarity
    
        
    
    import mp.mp_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z,w,u) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))),u);


    %% test 1
    len1 = 500; r = 300; k= 2^14;

    bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
    bar2 = [imgaussfilt(normrnd(0,1,len2,1),2.3); bar1; imgaussfilt(normrnd(0,1,len2,1),2.3)];

    bit1 = ones(len1,1);
    
    if nargin >=2
        bit1(1:bitmask) = 0;
        bit1(end-bitmask+1:end) = 0;
    end
    
    tic
    [mp, mpI,mpD] = comparisonFun(bar1,bar2,bit1,r,islinear);
    toc
    
    len2Full = length(bar2);
    
    % pos should be pos=1001
    import CBT.Hca.UI.Helper.get_best_parameters_mp;
    [ maxcoef,pos,or,idxpos] = get_best_parameters_mp( mp,mpI,mpD, 1, 50,islinear,len2Full);

% %     
% pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
% 
% if ~or(1)
%     frag1 = flipud(bar1);
% else
%     frag1 = bar1;
% end
% % 
% frag2 = bar2(pos(1):pos(1)+length(bar1)-1);
% 
% bit1 = find(bit1,1,'first');
% 
% frag1subseq = frag1(idxpos-bit1+1:idxpos-bit1+1+r-1);
% frag2subseq = frag2(idxpos-bit1+1:idxpos-bit1+1+r-1);
% 
% pcc(frag1,frag2)
% %     figure,plot(bar2)
%     hold on
%     plot(pos(1):pos(1)+length(bar1)-1,bar1)
%     frag1 = bar1;
%     frag2 = bar2(pos(1):pos(1)+length(bar1)-1);
%     pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
%     pcc(frag1,frag2)

%     assert(pos(1) == len2+1)


%     assert(pos(1) == len2+1)

end

