function [] = plot_best_bar(fig1, barcodeGen, consensusStruct, comparisonStruct,theoryStruct, maxcoef,userDefinedSeqCushion )
    % plot_best_bar
    
    % plots best barcode vs theory in case barcode is always larger than
    % theory
    if nargin < 7
        userDefinedSeqCushion = 0;
    end

    len1=length(barcodeGen);
    % barcode orientations
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');

    % number and value of the best barcode
    [dd,ii] =max(maxcoef(:,1));
    
    % load theory file
    fileID = fopen(theoryStruct{comparisonStruct{ii}.idx}.filename,'r');
    formatSpec = '%f';
    theorBar = fscanf(fileID,formatSpec);
    fclose(fileID);
    
    niceName = theoryStruct{comparisonStruct{ii}.idx}.name;
    pl = [strfind(niceName,'NC') strfind(niceName,'NZ')];
    niceName = niceName(pl:end);
    pl = [strfind(niceName,'|') strfind(niceName,' ')];
    niceName = strrep(niceName(1:(min(pl)-1)),'_','\_');
    if isempty(niceName)
        niceName = strrep(theoryStruct{comparisonStruct{ii}.idx}.name,'_','\_');
    end
        
    
    % theory length
    thrLen = theoryStruct{comparisonStruct{ii}.idx}.length;
    
    % bitmask. In case of linear barcode, would like to modify this
    theorBit = ones(1,thrLen);
    
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
    
    
    % expBar with expanded cushion
    expBar = [ repmat(nan,1,userDefinedSeqCushion) expBar repmat(nan,1,userDefinedSeqCushion)];

%     numSt = min(userDefinedSeqCushion,userDefinedSeqCushion-comparisonStruct{ii}.pos(1));
    % don't allow theory start to loop over
    theoryStart = comparisonStruct{ii}.pos(1)-userDefinedSeqCushion;
    theoryEnd = comparisonStruct{ii}.pos(1)-userDefinedSeqCushion+length(expBar)-1;
    
    if theoryStart < 1
        % in circular case, these should be taken from the end of theorBar
        theorBar = [ repmat(nan,abs(theoryStart)+1,1); theorBar];
        theoryEnd = theoryEnd + abs(theoryStart)+1;
        theoryStart = 1;
        thrLen = thrLen+abs(theoryStart)+1;
    end
    
    if theoryEnd > thrLen % now split this into linear and nonlinear case..
        theorBar = [ theorBar; theorBar(1:theoryEnd-thrLen)];
    end
    

	barStruct.bar1 = expBar;
    barStruct.bar2= theorBar;

    if comparisonStruct{ii}.or(1) == 1
        barStruct.matchTable = [1 length(expBar) theoryStart theoryEnd ...
            comparisonStruct{ii}.or(1)  ];
    else
        barStruct.matchTable = [1 length(expBar) theoryEnd theoryStart ...
           comparisonStruct{ii}.or(1)  ];
    end
    
    import CBT.Hca.UI.Helper.create_full_table;
    [temp_table,barfragq, barfragr] = create_full_table(barStruct.matchTable, barStruct.bar1,barStruct.bar2,1);

    barfragq{1}(~isnan(barfragq{1})) = zscore(barfragq{1}(~isnan(barfragq{1})));
    barfragr{1}(~isnan(barfragr{1})) = zscore(barfragr{1}(~isnan(barfragr{1})));

%     figure,
    plot(barfragq{1})
    hold on
    plot(barfragr{1})
    
  

    xlabel('Position along the sequence cushion (px)','Interpreter','latex')
    ylabel('Z-scored','Interpreter','latex')
    if ii <= len1
        name = num2str(ii);
    else
        name = 'consensus';
    end

    title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
    %
    legend({strcat(['$\hat C_{\rm ' name '}=$' num2str(dd,'%0.2f')]), niceName},'Interpreter','latex')



%     
%     
%     % flip barcode if it gabe best match when flipped                  
%     if orientation(ii,1) == 2
%         expBar = fliplr(expBar);
%         expBit = fliplr(expBit);
%     end
%     % fit positions
%     fitPositions = pos(ii,1):pos(ii,1)+length(expBar)-1;
%     
%     if min(fitPositions) <= 0 || max(fitPositions) > thrLen
%         %<= 0 should not happen
% %         max = 
% %         plot(repmat(length(theorBar),1,length(-3:3)),-3:3, 'black');
% %         hold on
%         % the match at best position loops around circularly. Shift the
%         theorBar = [theorBar; theorBar(1:end-1) ];
%         theorBit = [theorBit theorBit(1:end-1) ];
% %         warning('Experiment vs. theory not printed');
%     end
%     
%     % position on the theory barcode
%     barFit = theorBar(fitPositions);
%     barBit = theorBit(fitPositions);
% 
% 
%     % mean and std of theory where bitmask of experiment is nonzero.
%     % if there were some zero's in the theory, then this would not be
%     % correct?
%     m1 = mean(barFit(logical(expBit)));
%     s1= std(barFit(logical(expBit)));
% 
%     % mean and std of experiment
%     m2 = mean(expBar(logical(expBit)));
%     s2= std(expBar(logical(expBit)));
% 
%     % check if score is the same
%     % depends if we have zscore 1 or zscore 0
% %     zscore(barFit(logical(expBit)),1)'*zscore(expBar(logical(expBit)),1)'/length(expBar(logical(expBit)))
% %     dd
% 
%     plot(fitPositions, ((expBar-m2)/s2) *s1+m1)
%     hold on
%     plot(fitPositions, barFit)
%     xlim([min(fitPositions) max(fitPositions)])
% 
%     xlabel('Position (px)','Interpreter','latex')
%     ylabel('Rescaled to theoretical intesity','Interpreter','latex')
%     if ii <= len1
%         name = num2str(ii);
%     else
%         name = 'consensus';
%     end
% 
%     title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%      %
%     legend({strcat(['$\hat C_{\rm ' name '}=$' num2str(dd,'%0.2f')]), niceName},'Interpreter','latex')

end

