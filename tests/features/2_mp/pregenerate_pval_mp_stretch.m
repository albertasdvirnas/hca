function [cthreshFun,params,sol,leng] = pregenerate_pval_mp_stretch(low,high,w,numBars,sets)


    
%% test which is plotted in paper, generating p-value table for len1 vs len2

sets.svList = zeros(1,numBars);
% generate sample linear structural variations
import Rand.generate_linear_sv;
% [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(low, high, x,sets),sets.svList,'UniformOutput',false);

% import mp.mp_profile_stomp_dna;
kk = 2^16;

import mp.mp_full_masked_profile_stomp_dna;
comparisonFun = @(x,y,z1,z2,w,u) mp_full_masked_profile_stomp_dna(x,y,z1,z2,w,kk,u);

% [mpAB,mpIAB,mpDAB] = comparisonFun(X1, X2, bitX1, bitX2, w,islinear);
    
    
leng = low-5:10:high+5;
% leng2 = low:20:high;
stretch = sets.stretch;
import Rand.generate_linear_sv;

sol = zeros(length(leng),length(leng));
params = cell(length(leng),length(leng));
for j = 1:length(leng)
    j
    for k=j:length(leng)
        [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(leng(j), leng(k), x,sets),sets.svList,'UniformOutput',false);
        ccMax = zeros(1,length(bar1));

%             ccMax2 = zeros(1,length(bar1));
            for idx =1:length(bar1)
                maxL = zeros(1,length(stretch));
                for l=1:length(stretch)
                    bar1stretch = imresize(bar1{idx},[1 round(leng(j)*stretch(l))]);
                    % stretch 
                    % run this for lengthJ

                     [mp,mpI]=  comparisonFun(bar1stretch', bar2{idx}', ones(length(bar1stretch),1), ones(length(bar2{idx}),1), w,~sets.circ);
%                   
                    maxL(l) = max(mp);
                end
                ccMax(idx) = max(maxL);
    
            end


        import Pvalue.compute_evd_params; % x0 should be better bound, make an analysis
        params{j,k} = compute_evd_params(ccMax(:),100);
        p = @(x) 1-(0.5+0.5*(1-betainc((x).^2,0.5,params{j,k}(1)/2-1,'upper'))).^params{j,k}(2) ;
        % we want to have approx region
        sol(j,k) = fzero(@(x) p(x)-0.01, [0 1]);
    end

end

sol2 = sol;

sol2 = sol2 + transpose(triu(sol2))-tril(sol2);

cthreshFun = @(x,y) interp2(leng, leng, sol2, x,y);
% cthreshFun = @(x,y) interp2(leng, leng, sol, x,y,'linear','extrap');


end

