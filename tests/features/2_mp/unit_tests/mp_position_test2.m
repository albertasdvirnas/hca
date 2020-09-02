function [pos] = mp_position_test2(len2,bitmask, len1,islinear)
    % run_mp_generation_tests
    % MP generation tests to see if Matrix Profile method runs properly (from both computational 
    % and results perspective ) for
    % sequence similarity
    
        
    if nargin < 3
        len1 = 500;
    end
    
    import mp.mp_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z,w,u) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))),u);


    %% test 1
    %len1 = 500;
    r = 300; k= 2^14;

    bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
    bar2 = [bar1(101:len1); imgaussfilt(normrnd(0,1,len2,1),2.3); bar1(1:100)];

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
    [ maxcoef,pos,or ] = get_best_parameters_mp( mp,mpI, mpD, 1, 50, islinear, len2Full);

%     assert(pos(1) == len2+1)

end

