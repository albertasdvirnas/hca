function [ sets ] = get_theory_sets( sets )
    % get_theory_sets
    %
    % This function asks for all the required settings for the generation
    % of the theory at this point.
    %
	%     Args:
    %         sets (struct): Input settings to the method
    % 
    %     Returns:
    %         sets: Return structure
    
    prompt = {'Choose nm/bps','concNetropsin Mmolar',...
        'concYOYO1 Mmolar','conc DNA Mmolar', 'psfSigmaWidth nm',...
        'pixelWidth_nm','deltaCut','isLinearTF','widthSigmasFromMean',...
        'compute free concentrations','model'};
    
    title = 'Parameter choice for barcode generation';
    dims = [1 35];
    definput = {'0.3','6','0.02','0.2','300','130','3','0','4','1','literature'};
    answer = inputdlg(prompt,title,dims,definput);

    
    if ~isempty(answer)
        sets.meanBpExt_nm  = str2double(answer{1});
        sets.concN=str2double(answer{2});
        sets.concY=str2double(answer{3});
        sets.concDNA= str2double(answer{4});
        sets.psfSigmaWidth_nm= str2double(answer{5});
        sets.pixelWidth_nm = str2double(answer{6});
        sets.deltaCut = str2double(answer{7});
        sets.isLinearTF = str2double(answer{8});
        sets.widthSigmasFromMean = str2double(answer{9});
        sets.computeFreeConcentrations = str2double(answer{10});
        sets.model = answer{11};
    else
        disp('Default theory settings are being used');
    end

end

