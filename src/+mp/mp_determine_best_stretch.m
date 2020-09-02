function [bestStretch,maxCoef] = mp_determine_best_stretch(query_txt,data_txt,sets)
    %
    if nargin < 3
        sets.stretch =  0.7:0.01:1.3;
        sets.r = 100;
        sets.k = 2^14;
        
    end
    
    stretch = sets.stretch;
    
    import mp.mp_masked_profile_stomp_dna;
    comparisonFun = @(x,y,z,w,u) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))),u);

    %% test 1
     r = sets.r; 
     k = sets.k;

     lenD = length(query_txt);
     query_bit = ones(1,lenD);
     len2 = length(data_txt);

    % possible stretching factors
%     stretch = 0.7:0.01:1.3;
%     
    % original interval
    t = 1:lenD;
%     import CBT.Hca.UI.Helper.get_best_parameters_mp;


    % pos should be pos=1001
    import CBT.Hca.UI.Helper.get_best_parameters_mp;
    import Comparison.interpolate_data;

    comparisonStruct = cell(1,length(stretch));
    
    % make a loop over stretch factors
     for j=1:length(stretch)
         % rescale
        len1 = round(lenD*stretch(j));
        y = interpolate_data(query_txt,lenD,lenD*stretch(j),sets);

%         y = imresize(query_txt, [1,len1]);
        % rescale bitmask (though for this it is assumed bitmask always
        % ones)
        query_bit(query_bit==0)=nan;
        bit1res = interp1(query_bit, linspace(1,lenD,lenD*stretch(j)));
        bit1res(isnan(bit1res))=0;

        if length(y) >= r
            [mp, mpI,mpD] = comparisonFun(y',data_txt',bit1res',r,~sets.circ);
    %         [mpExample,mpIExample] = mp_profile_stomp_dna(y', bTS',r,kk);
        else
            mp = [];
            mpI = [];
        end

        [ comparisonStruct{j}.maxcoef,comparisonStruct{j}.pos,comparisonStruct{j}.or ] =...
            get_best_parameters_mp( mp, mpI, mpD, 1, 50, ~sets.circ, len2);

        % save results
        comparisonStruct{j}.mp = mp;
        comparisonStruct{j}.mpI = mpI;

     end

     %% max coeff based on stretch factor
    maxCoef = cellfun(@(x) x.maxcoef,comparisonStruct);
%     figure,plot(stretch,maxCoef)

    [a, b]  = max(maxCoef);
    bestStretch = stretch(b);
%     sets.analysis.stretchFactors = bestStretch-0.05:0.01:bestStretch+0.05;
    % 
    
end

