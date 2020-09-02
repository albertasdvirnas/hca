function [dataCC,t] = precompute_pvalue_files_multi_theory(setsFile,lengths )
    % precompute_pvalue_files_multi_theory
    
    % we allow there to be multi theories of different lengths.
    % this way it is better to include case when barcodes are circular.
    
    % Pregenerates p-value database file
    %     Args:
    %         fullPath, params 
    % 
    %     Returns:
    % 
    %     Example:
    %
    
    params = ini2struct(setsFile);
 
        
    tic
    % compute the psf
    psf = params.psfSigmaWidth_nm/params.pixelWidth_nm;
    
    % here allow to choose a method?
    

    % the short lengths that we need to compute for
    comparisonFun = @(x,y,z,w) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));

    import CBT.Hca.Core.Pvalue.compute_random_max_cc;
%     import SignalRegistration.unmasked_pcc_corr;
	import CBT.Hca.Import.load_pval_struct;
    import CBT.Hca.Core.Pvalue.convolve_bar;

    % lengths to compute
%     lens =params.lenMin:params.lenMax;
    
    lenMin=params.lenMin;
    lenMax=params.lenMax;
    isLinearQ = params.islinearQ;
    isLinearD = params.islinearD;

%     lengths = params.lengths;
    % 
%     vals = zeros(1,params.lenMax);
%     data = cell(1,params.lenMax);
    
    dataCC = cell(1,length(lengths));
    
    numRnd = params.numRnd;
    for k=1:length(lengths)
        data = cell(1,lenMax);
%         k

        % compute the long random barcode
        if isLinearD
            rand2 = normrnd(0,1, 1, 2*lengths(k));
            rand2 = convolve_bar(rand2, psf, 2*lengths(k) );
            rand2 = rand2(1:lengths(k));
        else %circular
            rand2 = normrnd(0,1, 1, lengths(k));
            rand2 = convolve_bar(rand2, psf, lengths(k) );
        end
        
        parfor i = lenMin:lenMax
            bit1 = ones(1,i);
            % we load the data every time, maybe change this with simpler
            % approach where we save the data only in the end

            % compute numRnd random barcodes for length i. 
            % note that in case it's linear, it would suffice to compute once
            for j = 1:numRnd
                % generate random barcode, cases: linear/circular
                if isLinearQ
                    rand1 = normrnd(0,1, 1, 2*i);
                    rand1 = convolve_bar(rand1, psf, 2*i );
                    rand1 = rand1(1:i);
                else %circular
                    rand1 = normrnd(0,1, 1, i);
                    rand1 = convolve_bar(rand1, psf, i );
                end

                % also should include stretching here, but that adds up another
                % loop
                ccM = comparisonFun(rand1, rand2, bit1);
%                 [ccM] = unmasked_pcc_corr(rand1, rand2, bit1);
                data{i}(j) = max(ccM(:));
            end    

        end
        dataCC{k} = data;
    end
    
    t=toc;
        %             [ccM] = unmasked_pcc_corr(rand1, rand2, bit1);
        %             data(i) = max(ccM(:));
        
%         disp(strcat(['Computing p-value for barcodes of length ' num2str(lenCur) ', already done ' num2str(i-1) ' out of ' num2str(length(lens))]));
%         data{i} = compute_random_max_cc(i,rand2, psf, numRnd);%%    

%         vals(end+1) = lenCur;
%         data{end+1} = dataNew';
% %         % checks if this length was already considered before         

%     
%     [ ~, data2 ] = load_pval_struct(fullPath);
%     vals2 = 1:max(length(data),length(data2));
%     for i =1:max(length(data),length(data2))
% %         [a,b] = ismember(lenCur, vals);
%         if ~isempty(data{i})
%             data2{i} =[data2{i};data{i}];
%         end
%         
%     end
% 
%         
%     % exports structure
%     import CBT.Hca.Export.export_pval_struct;
%     export_pval_struct( fullPath,vals2,data2 );
        
end

