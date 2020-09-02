function [cthreshFun,params,sol,leng] = pregen_pval_mp(low,high,w,numBars,sets)


    
%% test which is plotted in paper, generating p-value table for len1 vs len2

sets.svList = zeros(1,numBars);
% generate sample linear structural variations
import Rand.generate_linear_sv;
% [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(low, high, x,sets),sets.svList,'UniformOutput',false);

import mp.mp_profile_stomp_dna;
leng = low-5:50:high+5;
% leng2 = low:20:high;

import Rand.generate_linear_sv;

sol = zeros(length(leng),length(leng));
params = cell(length(leng),length(leng));
for j = 1:length(leng)
    j
    for k=j:length(leng)
        [bar1,bar2,matchTable,lengths]  = arrayfun(@(x) generate_linear_sv(leng(j), leng(k), x,sets),sets.svList,'UniformOutput',false);

        ccMax = zeros(1,length(bar1));
        for idx =1:length(bar1)

            % run this for lengthJ
            if sets.circ
                [mp,mpI] = mp_profile_stomp_dna([bar1{idx} bar1{idx}(1:leng(j)-1)]', bar2{idx}',w,2^16);
            else
                [mp,mpI] = mp_profile_stomp_dna([bar1{idx}]', bar2{idx}',w,2^16);
            end 

            ccMax(idx) = max(mp);

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

