function [ p, cal ] = compute_p_val_score(cMaxVals, pvalData, barLen,strFac )
    % compute_p_val_score
    %
    % Compute actual score
    %     Args:
    %         cMaxVals, pvalData, barLen,strFac
    % 
    %     Returns:
    %         p, cal
    %     Example:
    %
    % maximum p-values

    p = zeros(1,length(barLen));
    cal = ones(1,length(barLen));
    for i=1:length(barLen)
        disp(strcat(['Computing p-values for barcode nr.' num2str(i)])); 
        % compute barcode lengths when stretched
        indx = round(barLen(i)*strFac);
        % intersect with corresponding values from the database
        [~, idxIntoA] = intersect(pvalData.len1, indx);
        if ~isequal(length(idxIntoA),length(indx))
            warning('Not enough data in the p-value database for this barcode');  
            p(i) = 1;
            cal(i) = 0;
        else
            % extract the cc coefficients for all the different
            % possible lengths
            
            if min(size(cell2mat(pvalData.data(idxIntoA)'))) > 1
                maxCCVals = max(cell2mat(pvalData.data(idxIntoA)'));
            else
                maxCCVals = cell2mat(pvalData.data(idxIntoA)');
            end
            
            import CBT.Hca.Core.Pvalue.compute_evd_params;
            params = compute_evd_params(maxCCVals(:),barLen(i)/4);

            p(i) = 1-(0.5+vpa(0.5)*(vpa(1,16)-vpa(betainc(cMaxVals(i).^2,0.5,params(1)/2-1,'upper'),16))).^params(2);
        end
    end
    

% 
% % also plot distribution
% 
% cc = linspace(min(maxCCVals),max(maxCCVals),1001);
% 
% f = figure;
% M = floor(sqrt(length(maxCCVals)));
% histogram(maxCCVals,'Normalization','pdf','NumBins',M);
% 
% hold on
% import Pvalue.functional_pdf;
% exactpdf = functional_pdf(cc,params);
% plot(cc,exactpdf,'black', 'linewidth',3)
% ylabel('PDF')
% xlabel('CC')
% title('Parametric form pdf fit for # num samples 1000')
% legend({'Correlation coefficient hist', 'Parametric PDF'})
% % saveas(f,fullfile('resultData','dist.eps'),'epsc')
    
end

