function [] = plot_best_concentric_image(hAxis, barcodeGen, consensusStruct, comparisonStruct,theoryStruct, maxcoef ,sets)
    % plot_best_concentric_image
 
    % plots best barcode vs theory as concentric plot

    % barcode orientations
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');
    
    % number and value of the best barcode. Might want to plot not the
    % first feature, but second feature
    [CCMAX,ii] =max(maxcoef(:,1));
    
    % load theory file
    fileID = fopen(theoryStruct{comparisonStruct{ii}.idx}.filename,'r');
    formatSpec = '%f';
    theorBar = fscanf(fileID,formatSpec)';
    fclose(fileID);
    sets.theory.isLinearTF = theoryStruct{comparisonStruct{ii}.idx}.isLinearTF;

    % theory length
    thrLen = theoryStruct{comparisonStruct{ii}.idx}.length;
    

    % load either theory barcode or the consensus barcode
    try
        expBar = barcodeGen{ii}.rawBarcode;
        expBit = barcodeGen{ii}.rawBitmask;
    catch
        try
            expBar = consensusStruct.rawBarcode;
            expBit = consensusStruct.rawBitmask;  
        catch
            expBar = consensusStruct{ii-length(barcodeGen)}.rawBarcode;
            expBit = consensusStruct{ii-length(barcodeGen)}.rawBitmask;         
        end
    end
    expLen = length(expBar);

    % interpolate to the length which gave best CC value
    expBar = interp1(expBar, linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch));
    expBit = expBit(round(linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch)));
    expBar(~expBit)= nan;
    expBar(expBit) = zscore(expBar(expBit));
    barStruct.bar1 = expBar;

    
    
    
%     
%     
%     niceName = theoryStruct{comparisonStruct{ii}.idx}.name;
%     pl = [strfind(niceName,'NC') strfind(niceName,'NZ')];
%     niceName = niceName(pl:end);
%     pl = [strfind(niceName,'|') strfind(niceName,' ')];
%     niceName = strrep(niceName(1:(min(pl)-1)),'_','\_');
%     if isempty(niceName)
%         niceName = theoryStruct{comparisonStruct{idx}.idx}.name;
%     end
        
    
    barStruct.bar2 = zscore(theorBar,1);

    if comparisonStruct{ii}.or(1) == 1
        barStruct.matchTable = [1 length(barStruct.bar1) comparisonStruct{ii}.pos(1) comparisonStruct{ii}.pos(1)+length(barStruct.bar1)-1 ...
           comparisonStruct{ii}.or(1)  ];
    else
        barStruct.matchTable = [1 length(barStruct.bar1) comparisonStruct{ii}.pos(1)+length(barStruct.bar1)-1 comparisonStruct{ii}.pos(1) ...
           comparisonStruct{ii}.or(1)  ];
    end
    

%     
sets.A = 'b';
sets.B = 'b';
    import CBT.Hca.UI.Helper.plot_concetric;
    plot_concetric(hAxis,barStruct,sets);
%         plot_concetric(hAxis,barStruct,sets);

%    concentric_plot(barStruct)


   
%    
%     
%     % flip barcode if it gives best match when flipped                  
%     if orientation(ii,1) == 2
%         expBar = fliplr(expBar);
%         expBit = fliplr(expBit);
%     end
%     % fit positions
%     x = pos(ii,1):pos(ii,1)+length(expBar)-1;
%     
%     if min(x) > 0 && max(x) < thrLen
%         % position on the theory barcode
%         thr = theorBar(x);
% 
%         % mean and std of theory where bitmask of experiment is nonzero.
%         % if there were some zero's in the theory, then this would not be
%         % correct?
%         m1 = mean(thr(logical(expBit)));
%         s1= std(thr(logical(expBit)));
% 
%         % mean and std of experiment
%         m2 = mean(expBar(logical(expBit)));
%         s2= std(expBar(logical(expBit)));
% 
%         bar = ((expBar-m2)/s2) *s1+m1;
%         
%         imshow([repmat(bar,50,1); repmat(thr,50,1)],[])
%         
% %         plot(fitPositions, ((expBar-m2)/s2) *s1+m1)
% %         hold on
% %         plot(fitPositions, barFit)
% %         xlim([min(fitPositions) max(fitPositions)])
% 
% %         xlabel('Position (px)','Interpreter','latex')
% %         ylabel('Rescaled to theoretical intesity','Interpreter','latex')
% %         if ii <= len1
% %             name = num2str(ii);
% %         else
% %             name = 'consensus';
% %         end
% % 
% %         title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
% %         legend({strcat(['$\hat C_{\rm ' name '}=$' num2str(dd,'%0.2f')]), matlab.lang.makeValidName(theoryStruct{comparisonStruct{ii}.idx}.name)},'Interpreter','latex')
%     else
%         warning('Experiment vs. theory not printed');
%     end
end

